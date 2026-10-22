/*
컴프2 중간과제 대체 팀플 2조 코드

조장: 아이디어, 기획, 일정, pm, 데이터 수집, 코드 검증 및 정리, 공유 폴더 관리, 코딩(파일 입출력, R파 검출, 메뉴, 터닝포인트)
팀원 A: 코딩(다운샘플링, 양자화)
팀원 B: 코딩(RLE, 델타 부호화)
팀원 C: 발표, 데이터 시각화(손실)
팀원 D: 발표, 데이터 시각화(무손실)

파일 작동법: 

 * [메뉴]
 *   1  MIT-BIH 기록 읽기 (기록 번호 입력: 100, 103, 119, 203, 208 있음, 골라서 입력해주세요)
 *   2  BIOPAC 측정 파일 읽기 (구간 선택, 데이터 중간에 심한 노이즈가 있어서 제거하고 2개로 쪼갬)
 *   3  압축 실험 → R파 검출 정확도(MIT-BIH만) + 결과표 + 방법별 최대 압축 + results_데이터이름.csv 저장
 *   0  종료
 *   사용 순서 예: 1 → 100 → 3 → 0   (3번은 1·2번으로 데이터를 먼저 읽어야 함)
 
 !! 데이터 파일은 data/ 폴더에 두고, 레포 최상위 폴더에서 실행해주세요 (Visual Studio는 프로젝트 폴더에 data/ 폴더를 두기) !!
 */


//공통 정의
#ifndef ECG_H
#define ECG_H

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_N      650000
#define MAX_BEATS  10000
#define MAX_RESULT 100

#define MITBIH_DIR   "data/"             //MIT-BIH 데이터 폴더 (실행 위치 기준)
#define BIOPAC_FILE  "data/ecg_t3_03.txt" //BIOPAC 측정 파일 (연구실 데이터라 레포에는 없음)

extern short orig[MAX_N];
extern short rec[MAX_N];
extern int   n_samples;
extern int   fs;

extern short         rle_val[MAX_N];
extern unsigned char rle_cnt[MAX_N];
extern int           n_runs;
extern char          diff[MAX_N];
extern short         kept[MAX_N];
extern int           n_kept;

extern int truth[MAX_BEATS];
extern int n_truth;
extern int beats_orig[MAX_BEATS];
extern int n_beats_orig;
extern int beats_rec[MAX_BEATS];
extern int n_beats_rec;

typedef struct {
    char   name[64];
    long   orig_bytes;
    long   comp_bytes;
    double cr;
    double prd;
    int    tp, fn, fp;
    double se, ppv;
    double hr_orig, hr_rec;
    int    success;
} Result;

extern Result results[MAX_RESULT];
extern int    n_results;

extern char data_name[64];
extern char data_tag[64];

int  read_mitbih(const char path[]);
int  read_truth(const char path[]);
int  read_biopac(const char path[], int col, int rate,
                 double start_sec, double end_sec);
void save_results_csv(const char path[]);

long rle_save(void);
void rle_restore(void);
long diff_save(void);
void diff_restore(void);

void find_min_max(const short x[], int n, int *pmin, int *pmax);
long skip_save(int k);
void skip_restore(int k);
long quant_save(int s);
void quant_restore(int s);

long tp_save(int levels);
void tp_restore(int levels);

int find_beats(const short x[], int n, int rate, int beats[]);
int count_match(int ref[], int nref, int test[], int ntest, int tol);

double calc_prd(void);
double calc_hr(int beats[], int nb, int rate);
void   evaluate(const char name[], long comp_bytes);
void   print_results(void);

#endif

//파일 입출력: mit-bih, biopac 읽고, 결과 csv 저장
int read_mitbih(const char path[])
{
    FILE *fp;
    char  line[256];
    int   col = -1;
    int   comma = 0;
    int   vals[3];
    int   i, cnt;

    fp = fopen(path, "r");
    if (fp == NULL)
        return -1;

    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }
    for (i = 0; line[i] != '\0'; i++) {
        if (line[i] == ',')
            comma++;
        if (line[i] == 'M' && line[i + 1] == 'L' && line[i + 2] == 'I' && line[i + 3] == 'I') {
            col = comma;
            break;
        }
    }
    if (col < 1 || col > 2) {
        printf("이 기록에는 MLII 열이 없습니다 (102, 104번은 사용하지 않음)\n");
        fclose(fp);
        return -1;
    }

    n_samples = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        cnt = sscanf(line, "%d,%d,%d", &vals[0], &vals[1], &vals[2]);
        if (cnt <= col)
            continue;
        if (n_samples >= MAX_N) {
            printf("샘플이 %d개를 넘어서 그만 읽습니다\n", MAX_N);
            break;
        }
        orig[n_samples] = (short)(vals[col] - 1024);
        n_samples++;
    }
    fclose(fp);

    if (n_samples == 0)
        return -1;
    fs = 360;
    return n_samples;
}

int is_beat_symbol(char c)
{
    char beat_symbols[] = "NLRBAaJSVrFejnE/fQ?";
    int  i;

    for (i = 0; beat_symbols[i] != '\0'; i++)
        if (c == beat_symbols[i])
            return 1;
    return 0;
}

int read_truth(const char path[])
{
    FILE *fp;
    char  line[256];
    char  time[32], sym[16];
    int   pos;

    fp = fopen(path, "r");
    if (fp == NULL)
        return -1;

    fgets(line, sizeof(line), fp);

    n_truth = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (sscanf(line, "%31s %d %15s", time, &pos, sym) != 3)
            continue;
        if (sym[1] == '\0' && is_beat_symbol(sym[0])) {
            if (n_truth >= MAX_BEATS)
                break;
            truth[n_truth] = pos;
            n_truth++;
        }
    }
    fclose(fp);
    return n_truth;
}

int read_biopac(const char path[], int col, int rate, double start_sec, double end_sec)
{
    FILE  *fp;
    char   line[512];
    double v[6];
    long   k = 0;
    long   start_k, end_k;
    int    cnt;

    if (col < 1 || col > 5)
        return -1;

    fp = fopen(path, "r");
    if (fp == NULL)
        return -1;

    start_k = lround(start_sec * rate);
    end_k   = (end_sec > 0) ? lround(end_sec * rate) : -1;

    n_samples = 0;
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (!((line[0] >= '0' && line[0] <= '9') || line[0] == '-'))
            continue;

        cnt = sscanf(line, "%lf %lf %lf %lf %lf %lf", &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]);
        if (cnt < col + 1)
            continue;

        if (end_k >= 0 && k >= end_k)
            break;
        if (k >= start_k) {
            if (n_samples >= MAX_N) {
                printf("샘플이 %d개를 넘어서 그만 읽습니다\n", MAX_N);
                break;
            }
            orig[n_samples] = (short)lround(v[col] * 1000.0);
            n_samples++;
        }
        k++;
    }
    fclose(fp);

    if (n_samples == 0)
        return -1;
    fs = rate;
    return n_samples;
}

void save_results_csv(const char path[])
{
    FILE *fp;
    int   i;

    fp = fopen(path, "w");
    if (fp == NULL) {
        printf("결과 파일을 만들 수 없습니다: %s\n", path);
        return;
    }

    fprintf(fp, "\xEF\xBB\xBF");
    fprintf(fp, "데이터,방법,원본크기(B),압축크기(B),CR,PRD(%%),Se(%%),PPV(%%),HR원본(bpm),HR복원(bpm),성공\n");
    for (i = 0; i < n_results; i++) {
        fprintf(fp, "%s,%s,%ld,%ld,%.3f,%.3f,%.2f,%.2f,%.2f,%.2f,%s\n",
                data_name, results[i].name, results[i].orig_bytes, results[i].comp_bytes,
                results[i].cr, results[i].prd, results[i].se, results[i].ppv,
                results[i].hr_orig, results[i].hr_rec,
                results[i].success ? "O" : "X");
    }
    fclose(fp);
    printf("결과를 %s에 저장했습니다 (%d줄)\n", path, n_results);
}

//R파 검출, 박동 비교
int find_beats(const short x[], int n, int rate, int beats[])
{
    int    block = 2 * rate;
    int    gap   = (int)(0.2 * rate);
    int    last  = -gap;
    int    nb    = 0;
    int    start, end, i;
    int    max_val;
    double sum, mean, th;

    for (start = 0; start < n; start += block) {
        end = start + block;
        if (end > n)
            end = n;

        sum = 0.0;
        max_val = x[start];
        for (i = start; i < end; i++) {
            sum += x[i];
            if (x[i] > max_val)
                max_val = x[i];
        }
        mean = sum / (end - start);
        th   = mean + 0.6 * (max_val - mean);

        for (i = start; i < end; i++) {
            if (i == 0 || i == n - 1)
                continue;
            if (x[i] > th && x[i] > x[i - 1] && x[i] >= x[i + 1] && i - last >= gap) {
                if (nb >= MAX_BEATS)
                    return nb;
                beats[nb] = i;
                nb++;
                last = i;
            }
        }
    }
    return nb;
}

int count_match(int ref[], int nref, int test[], int ntest, int tol)
{
    int i;
    int j  = 0;
    int tp = 0;

    for (i = 0; i < nref; i++) {
        while (j < ntest && test[j] < ref[i] - tol)
            j++;

        if (j < ntest && test[j] <= ref[i] + tol) {
            tp++;
            j++;
        }
    }
    return tp;
}

//무손실

//RLE
long rle_save(void)
{
    int i;

    if (n_samples <= 0)
        return -1;

    rle_val[0] = orig[0];
    rle_cnt[0] = 1;
    n_runs = 1;

    for (i = 1; i < n_samples; i++) {
        if (orig[i] == rle_val[n_runs - 1] && rle_cnt[n_runs - 1] < 255) {
            rle_cnt[n_runs - 1]++;
        }
        else {
            rle_val[n_runs] = orig[i];
            rle_cnt[n_runs] = 1;
            n_runs++;
        }
    }
    return (long)n_runs * 3;
}

void rle_restore(void)
{
    int i, j;
    int pos = 0;

    for (i = 0; i < n_runs; i++) {
        for (j = 0; j < rle_cnt[i] && pos < MAX_N; j++) {
            rec[pos] = rle_val[i];
            pos++;
        }
    }
}

//델타 부호화

long diff_save(void)
{
    int i, d;

    if (n_samples <= 0)
        return -1;

    diff[0] = 0;
    for (i = 1; i < n_samples; i++) {
        d = orig[i] - orig[i - 1];
        if (d < -128 || d > 127) {
            printf("  [델타 부호화] %d번째 샘플의 차이 %d가 1바이트 범위(-128~127)를 넘어 건너뜁니다\n", i, d);
            return -1;
        }
        diff[i] = (char)d;
    }
    return 2 + (long)(n_samples - 1);
}

void diff_restore(void)
{
    int i;

    if (n_samples <= 0)
        return;

    rec[0] = orig[0];
    for (i = 1; i < n_samples; i++)
        rec[i] = (short)(rec[i - 1] + diff[i]);
}

//손실

//다운샘플링
static int quant_min = 0;

void find_min_max(const short x[], int n, int *pmin, int *pmax)
{
    int i;

    *pmin = x[0];
    *pmax = x[0];
    for (i = 1; i < n; i++) {
        if (x[i] < *pmin)
            *pmin = x[i];
        if (x[i] > *pmax)
            *pmax = x[i];
    }
}

long skip_save(int k)
{
    int i;

    if (k < 1 || k > n_samples)
        return -1;
    n_kept = 0;
    for (i = 0; i < n_samples; i += k)
        kept[n_kept++] = orig[i];
    return (long)n_kept * 2;
}

void skip_restore(int k)
{
    int j, m, i, a, step;

    if (k < 1 || k > n_samples || n_kept < 1)
        return;
    for (j = 0; j < n_kept; j++) {
        a = kept[j];
        if (j + 1 < n_kept)
            step = kept[j + 1] - a;
        else if (j > 0)
            step = a - kept[j - 1];
        else
            step = 0;
        for (m = 0; m < k; m++) {
            i = j * k + m;
            if (i >= n_samples)
                break;
            rec[i] = (short)(a + step * m / k);
        }
    }
}

//양자화

long quant_save(int s)
{
    int i, mx, top, b;

    if (s < 1 || s > 15 || n_samples < 1)
        return -1;
    find_min_max(orig, n_samples, &quant_min, &mx);
    for (i = 0; i < n_samples; i++)
        kept[i] = (short)((orig[i] - quant_min) >> s);
    n_kept = n_samples;
    top = (mx - quant_min) >> s;
    b = 0;
    while ((1 << b) <= top)
        b++;
    if (b == 0)
        b = 1;
    return ((long)n_samples * b + 7) / 8 + 2;
}

void quant_restore(int s)
{
    int i;

    if (s < 1 || s > 15)
        return;
    for (i = 0; i < n_samples; i++)
        rec[i] = (short)((kept[i] << s) + (1 << (s - 1)) + quant_min);
}

//터닝포인트

long tp_save(int levels)
{
    int   i, s1, s2, lv, n;
    short x0, x1, x2;

    if (levels < 1 || n_samples < 3)
        return -1;

    for (i = 0; i < n_samples; i++)
        kept[i] = orig[i];
    n_kept = n_samples;

    for (lv = 0; lv < levels; lv++) {
        n = n_kept;
        if (n < 3)
            return -1;

        x0 = kept[0];
        n_kept = 1;

        for (i = 1; i + 1 < n; i += 2) {
            x1 = kept[i];
            x2 = kept[i + 1];

            if (x1 > x0)      s1 = 1;
            else if (x1 < x0) s1 = -1;
            else              s1 = 0;

            if (x2 > x1)      s2 = 1;
            else if (x2 < x1) s2 = -1;
            else              s2 = 0;

            if (s1 * s2 < 0)
                x0 = x1;
            else
                x0 = x2;

            kept[n_kept] = x0;
            n_kept++;
        }
    }
    return (long)n_kept * 2;
}

void tp_restore(int levels)
{
    int  j, m, pos, a, b;
    int  step = 1 << levels;

    for (j = 0; j + 1 < n_kept; j++) {
        a = kept[j];
        b = kept[j + 1];
        for (m = 0; m < step; m++) {
            pos = j * step + m;
            if (pos >= n_samples)
                break;
            rec[pos] = (short)((a * (step - m) + b * m) / step);
        }
    }
    for (pos = (n_kept - 1) * step; pos < n_samples; pos++)
        rec[pos] = kept[n_kept - 1];
}

//평가, 결과 출력

double calc_prd(void)
{
    double mean = 0.0, err = 0.0, var = 0.0, d;
    int    i;

    if (n_samples < 1)
        return 0.0;

    for (i = 0; i < n_samples; i++)
        mean += orig[i];
    mean = mean / n_samples;

    for (i = 0; i < n_samples; i++) {
        d = orig[i] - rec[i];
        err += d * d;
        d = orig[i] - mean;
        var += d * d;
    }
    if (var == 0.0)
        return 0.0;
    return 100.0 * sqrt(err / var);
}

double calc_hr(int beats[], int nb, int rate)
{
    double avg_gap_sec;

    if (nb < 2)
        return 0.0;
    avg_gap_sec = (double)(beats[nb - 1] - beats[0]) / (nb - 1) / rate;
    return 60.0 / avg_gap_sec;
}

void evaluate(const char name[], long comp_bytes)
{
    Result r;
    int    tp;

    if (n_results >= MAX_RESULT) {
        printf("결과표가 꽉 찼습니다 (%d줄)\n", MAX_RESULT);
        return;
    }

    n_beats_rec = find_beats(rec, n_samples, fs, beats_rec);
    tp = count_match(beats_orig, n_beats_orig, beats_rec, n_beats_rec, (int)(0.05 * fs));

    strncpy(r.name, name, sizeof(r.name) - 1);
    r.name[sizeof(r.name) - 1] = '\0';
    r.tp = tp;
    r.fn = n_beats_orig - tp;
    r.fp = n_beats_rec - tp;

    r.se  = (n_beats_orig > 0) ? 100.0 * tp / n_beats_orig : 0.0;
    r.ppv = (n_beats_rec  > 0) ? 100.0 * tp / n_beats_rec  : 0.0;

    r.prd     = calc_prd();
    r.hr_orig = calc_hr(beats_orig, n_beats_orig, fs);
    r.hr_rec  = calc_hr(beats_rec, n_beats_rec, fs);

    r.orig_bytes = (long)n_samples * 2;
    r.comp_bytes = comp_bytes;
    r.cr = (comp_bytes > 0) ? (double)r.orig_bytes / comp_bytes : 0.0;

    r.success = (r.se >= 99.0 && r.ppv >= 99.0 && fabs(r.hr_orig - r.hr_rec) < 1.0);

    results[n_results] = r;
    n_results++;
}

void print_padded(const char s[], int width)
{
    int i = 0, cols = 0;

    while (s[i] != '\0') {
        if ((unsigned char)s[i] >= 0xE0) {
            cols += 2;
            i += 3;
        } else {
            cols += 1;
            i += 1;
        }
    }
    printf("%s", s);
    while (cols < width) {
        printf(" ");
        cols++;
    }
}

void print_results(void)
{
    int i;

    printf("\n");
    print_padded("방법", 26);
    printf("    CR   PRD(%%)     Se(%%)   PPV(%%)   ΔHR(bpm)  성공\n");
    printf("------------------------------------------------------------------------------\n");
    for (i = 0; i < n_results; i++) {
        print_padded(results[i].name, 26);
        printf("%6.2f  %7.2f  %8.2f  %7.2f  %9.2f     %s\n",
               results[i].cr, results[i].prd, results[i].se, results[i].ppv,
               fabs(results[i].hr_orig - results[i].hr_rec),
               results[i].success ? "O" : "X");
    }
    printf("(원본 HR %.1f bpm, 성공 = Se·PPV 99%% 이상 그리고 ΔHR 1 bpm 미만)\n",
           n_results > 0 ? results[0].hr_orig : 0.0);
}

//전역변수처리

short orig[MAX_N];
short rec[MAX_N];
int   n_samples = 0;
int   fs = 360;

short         rle_val[MAX_N];
unsigned char rle_cnt[MAX_N];
int           n_runs = 0;
char          diff[MAX_N];
short         kept[MAX_N];
int           n_kept = 0;

int truth[MAX_BEATS];
int n_truth = 0;
int beats_orig[MAX_BEATS];
int n_beats_orig = 0;
int beats_rec[MAX_BEATS];
int n_beats_rec = 0;

Result results[MAX_RESULT];
int    n_results = 0;

char data_name[64] = "";
char data_tag[64]  = "";

int loaded = 0;

//입력처리

int read_int(void)
{
    int v, c, r;

    while ((r = scanf("%d", &v)) != 1) {
        if (r == EOF) return 0;
        while ((c = getchar()) != '\n' && c != EOF) ;
        printf("숫자를 입력하세요: ");
    }
    return v;
}

double read_double(void)
{
    double v;
    int c, r;

    while ((r = scanf("%lf", &v)) != 1) {
        if (r == EOF) return 0;
        while ((c = getchar()) != '\n' && c != EOF) ;
        printf("숫자를 입력하세요: ");
    }
    return v;
}

//메뉴

void menu_read_mitbih(void)
{
    int  num;
    char path[200];

    printf("기록 번호 (예: 100): ");
    num = read_int();

    sprintf(path, "%s%d.csv", MITBIH_DIR, num);
    if (read_mitbih(path) < 0) {
        printf("신호 파일을 읽지 못했습니다: %s\n", path);
        return;
    }

    sprintf(path, "%s%dannotations.txt", MITBIH_DIR, num);
    if (read_truth(path) < 0) {
        printf("정답 파일을 읽지 못했습니다: %s (정답 없이 진행)\n", path);
        n_truth = 0;
    }

    loaded = 1;
    sprintf(data_name, "MIT-BIH %d", num);
    sprintf(data_tag, "MITBIH_%d", num);
    printf("MIT-BIH %d: 샘플 %d개, %d Hz, 정답 박동 %d개\n", num, n_samples, fs, n_truth);
}

void menu_read_biopac(void)
{
    int    choice;
    double start_sec, end_sec;

    printf("구간 선택 (60.0~62.5초는 전원 잡음이라 제외)\n");
    printf("  1. 62.6 ~ 182.5초 (박동 약 146개, 추천)\n");
    printf("  2. 0 ~ 59.7초 (박동 약 72개)\n");
    printf("  3. 직접 입력\n");
    printf("선택: ");
    choice = read_int();

    if (choice == 1)      { start_sec = 62.6; end_sec = 182.5; }
    else if (choice == 2) { start_sec = 0.0;  end_sec = 59.7;  }
    else if (choice == 3) {
        printf("시작(초): ");
        start_sec = read_double();
        printf("끝(초, 0이면 끝까지): ");
        end_sec = read_double();
    }
    else { printf("1~3 중에서 고르세요.\n"); return; }

    if (read_biopac(BIOPAC_FILE, 2, 2000, start_sec, end_sec) < 0) {
        printf("파일을 읽지 못했습니다: %s\n", BIOPAC_FILE);
        return;
    }
    n_truth = 0;
    loaded = 1;
    if (end_sec > 0) {
        sprintf(data_name, "BIOPAC %.1f-%.1fs", start_sec, end_sec);
        sprintf(data_tag, "BIOPAC_%.1f-%.1fs", start_sec, end_sec);
    } else {
        sprintf(data_name, "BIOPAC %.1fs-end", start_sec);
        sprintf(data_tag, "BIOPAC_%.1fs-end", start_sec);
    }
    printf("BIOPAC 실험 3 유도 II: %.1f ~ %.1f초, 샘플 %d개, %d Hz\n",
           start_sec, end_sec, n_samples, fs);
}

void run_one(const char name[], long bytes, int method, int param)
{
    if (bytes < 0) {
        printf("  %-24s 건너뜀 (구현 전이거나 실패)\n", name);
        return;
    }
    switch (method) {
    case 1: rle_restore();        break;
    case 2: diff_restore();       break;
    case 3: skip_restore(param);  break;
    case 4: quant_restore(param); break;
    case 5: tp_restore(param);    break;
    }
    evaluate(name, bytes);
}

int last_success(int before)
{
    if (n_results == before)
        return 0;
    return results[n_results - 1].success;
}

void print_limit(const char label[], int from, int to)
{
    int i, best = -1;
    Result *r;

    for (i = from; i < to; i++) {
        if (!results[i].success)
            break;
        best = i;
    }

    if (best < 0)
        printf("  %s: 처음 설정부터 실패\n", label);
    else {
        printf("  ");
        print_padded(results[best].name, 24);
        printf("최대 → CR %.2f, PRD %.2f%%\n", results[best].cr, results[best].prd);
    }

    if (i < to) {
        r = &results[i];
        printf("    (멈춘 곳: %s —", r->name);
        if (r->se < 99.0)  printf(" Se %.2f%%", r->se);
        if (r->ppv < 99.0) printf(" PPV %.2f%%", r->ppv);
        if (fabs(r->hr_orig - r->hr_rec) >= 1.0)
            printf(" ΔHR %.2f bpm", fabs(r->hr_orig - r->hr_rec));
        printf(")\n");
    }
}

void menu_run_all(void)
{
    int  i, k, s, lv, tp, before, skip_from, skip_to, q_from, q_to, tp_from, tp_to;
    char name[64];
    char path[128];

    printf("데이터: %s\n", data_name);
    n_results = 0;
    n_beats_orig = find_beats(orig, n_samples, fs, beats_orig);
    printf("원본에서 찾은 박동: %d개\n", n_beats_orig);

    if (n_truth > 0 && n_beats_orig > 0) {
        tp = count_match(truth, n_truth, beats_orig, n_beats_orig, (int)(0.15 * fs));
        printf("R파 검출 정확도 (의사 판독 정답 %d개 대비): Se %.1f%%, PPV %.1f%%\n",
               n_truth, 100.0 * tp / n_truth, 100.0 * tp / n_beats_orig);
    }

    run_one("RLE", rle_save(), 1, 0);
    run_one("델타 부호화 (8비트)", diff_save(), 2, 0);

    skip_from = n_results;
    for (k = 2; k <= n_samples / 2 && n_results < MAX_RESULT - 20; k++) {
        sprintf(name, "다운샘플링 1/%d (%.0fHz)", k, (double)fs / k);
        before = n_results;
        run_one(name, skip_save(k), 3, k);
        if (!last_success(before))
            break;
    }
    skip_to = n_results;

    q_from = n_results;
    for (s = 1; s <= 15 && n_results < MAX_RESULT - 1; s++) {
        sprintf(name, "양자화 %d비트 제거", s);
        before = n_results;
        run_one(name, quant_save(s), 4, s);
        if (!last_success(before))
            break;
    }
    q_to = n_results;

    tp_from = n_results;
    for (lv = 1; lv <= 15 && n_results < MAX_RESULT - 1; lv++) {
        sprintf(name, "터닝포인트 %d단계 (CR %d)", lv, 1 << lv);
        before = n_results;
        run_one(name, tp_save(lv), 5, lv);
        if (!last_success(before))
            break;
    }
    tp_to = n_results;

    print_results();

    printf("\n최대 압축 (성공 = Se·PPV 99%% 이상, ΔHR 1 bpm 미만. PRD는 참고)\n");
    for (i = 0; i < n_results; i++)
        if (i < skip_from)
        {
            printf("  ");
            print_padded(results[i].name, 24);
            printf("%s → CR %.2f, PRD %.2f%%\n",
                   results[i].success ? "성공" : "실패", results[i].cr, results[i].prd);
        }
    print_limit("다운샘플링", skip_from, skip_to);
    print_limit("양자화", q_from, q_to);
    print_limit("터닝포인트", tp_from, tp_to);

    sprintf(path, "results_%s.csv", data_tag);
    save_results_csv(path);
}


int main(void)
{
    int choice;

    while (1) {
        printf("\n===== ECG 압축 분석기 =====\n");
        printf("1. MIT-BIH 기록 읽기\n");
        printf("2. BIOPAC 측정 파일 읽기\n");
        printf("3. 압축 실험 → 결과표 + CSV 저장\n");
        printf("0. 종료\n");
        printf("선택: ");
        choice = read_int();

        if (choice == 3 && (!loaded || n_samples == 0)) {
            printf("먼저 파일을 읽으세요 (1번 또는 2번)\n");
            continue;
        }

        switch (choice) {
        case 1: menu_read_mitbih();  break;
        case 2: menu_read_biopac();  break;
        case 3: menu_run_all();      break;
        case 0: printf("종료합니다.\n"); return 0;
        default: printf("0~3 중에서 고르세요.\n");
        }
    }
}
