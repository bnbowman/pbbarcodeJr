#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

#define M 128
#define N 128
#define MAX(x,y) (((x) > (y)) ? (x) : (y))
#define MAX3(a,b,c) (MAX(MAX(a,b), c))
#define MAX4(a,b,c,d) (MAX(MAX(a,b), MAX(c,d)))

int* allocate_dp_mat() {
    return (int*) calloc(N*M, sizeof(int));
}

int compute_align_score(int* dp_mat, char* tSeq, char* qSeq) {
    int ipenalty   = -7;
    int dpenalty   = -7;
    int match      =  4;
    int mpenalty   = -13;
    int bpenalty   = -4;
    int best_score = 0;
    int iscore     = 0;
    int dscore     = 0;
    int mscore     = 0;
    bool is_branch;
    int i,j;
    int max_i, max_j;

    memset(dp_mat, 0, M*N*sizeof(int));

    max_i = strlen(tSeq);
    max_j = strlen(qSeq);
    for (i = 1; i < strlen(tSeq) + 1; i++) {
	for (j = 1; j < strlen(qSeq) + 1; j++) {

	    iscore = dp_mat[i*M + j-1] + ipenalty;

            is_branch = (tSeq[j] == tSeq[i-1]) ? true : false;
	    dscore = dp_mat[(i-1)*M + j] + ((is_branch) ? bpenalty : dpenalty);

	    mscore = dp_mat[(i-1)*M + j-1] + ((tSeq[i-1] == qSeq[j-1]) ? match : mpenalty);

	    dp_mat[i*M + j] = MAX3(iscore, dscore, mscore);

 	    if (j == max_j && dp_mat[i*M + j] >= best_score) 
	        best_score = dp_mat[i*M + j];
	}
    }
    return best_score;
}

void compute_align_scores(int* scores, int n, int* dp_mat, char* tSeq, 
                          char** qSeqs) {
    int i = 0;
    for (i; i < n; i++) {
        scores[i] = compute_align_score(dp_mat, tSeq, qSeqs[i]);
    }
}


void print_dp_mat(int* dp_mat, char* tSeq, char* qSeq) {
    int i,j,score;
    for (i = 0; i < strlen(tSeq) + 1; i++) {
        if (i == 0) {
            printf("      %c  ", tSeq[0]);
        } else {
            printf("%c  ", tSeq[i]);
        }
    }
    printf("\n");

    for (j = 0; j < strlen(qSeq) + 1; j++) {
        if (j == 0) {
            printf("   ");
        } else {
            printf("%c  ", qSeq[j-1]);
        }

    	for (i = 0; i < strlen(tSeq) + 1; i++) {
    	    score = dp_mat[i*M + j];
    	    if (score >= 0 && score < 10) {
    	        printf("%d  ", score);
    	    } else if (score <= -10) {
    	        printf("%d", score);
    	    } else {
    	        printf("%d ", score);
            }
    	}
    	printf("\n");
    }
}
