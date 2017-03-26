#include "projet.h"
#include <mpi.h>

/* 2017-02-23 : version 1.0 */
#define FALSE 0
#define TRUE 1

unsigned long long int node_searched = 0;

struct recTree_t {
    int parentId;
    int move;
    tree_t *tree;
    result_t *result;
};
typedef struct recTree_t recTree_t;
void evaluate(tree_t * T, result_t *result);

void pre_evaluate (tree_t *T, result_t *result) {
    node_searched++;

    /* MPI vars */
    int nb_proc = 2;
    MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

    /* Tree handling
     * beg is the index of the beginning
     * of the current line of the tree */
    recTree_t *preEvalTrees;
    int beg, sizeTree=1;

    /* Array to store moves */
    int nb_tasks, n_moves;
    move_t moves[MAX_MOVES];

    /* Treee initial allocation,
     * it is then reallocated along the execution */
    if ( (preEvalTrees=malloc(sizeof(recTree_t))) == NULL) {
        fprintf(stderr,"malloc error in pre_evaluate()\n");
        exit(1);
    }
    /* Set structure up for first element */
    preEvalTrees[0].parentId = -1;
    preEvalTrees[0].tree = T;
    preEvalTrees[0].result = result;

    /* Original tree_t/result_t preparation phase
     * We skip tests for draw/victory, max depth
     * and n_moves because we assume real game is
     * given in input */
    result->score = -MAX_SCORE -1;
    result->pv_length = 0;
    compute_attack_squares(T);
    nb_tasks = n_moves = generate_legal_moves(T, moves);
    beg = sizeTree;
    sizeTree+=n_moves;
    if ( (preEvalTrees=realloc(preEvalTrees, sizeTree*sizeof(recTree_t))) == NULL) {
        fprintf(stderr, "realloc error in pre_evaluate()\n");
        exit(1);
    }
    // TODO some comments here?
    for (int i=beg ; i<sizeTree ; i++) {
        preEvalTrees[i].parentId = 0;
        preEvalTrees[i].move = moves[i-beg];
        preEvalTrees[i].tree = malloc(sizeof(tree_t));
        preEvalTrees[i].result = malloc(sizeof(result_t));
        if (!preEvalTrees[i].tree || !preEvalTrees[i].result) {
            fprintf(stderr,"malloc error in pre_evaluate()\n");
            exit(1);
        }
        play_move(preEvalTrees[0].tree, moves[i-beg], preEvalTrees[i].tree);
    }

    while (nb_tasks < 10*nb_proc) {
        int nbContinue=0, sum=0;
        int bound = beg+nb_tasks;
        for (int i=beg ; i<bound ; i++) {
            node_searched++;
            tree_t *tmpTree = preEvalTrees[i].tree;
            result_t *tmpResult = preEvalTrees[i].result;

            tmpResult->score = -MAX_SCORE-1;
            tmpResult->pv_length = 0;
            if (test_draw_or_victory(tmpTree, tmpResult)){
                nbContinue++;
                continue;
            }
            if (tmpTree->depth ==0) {
                tmpResult->score = (2*tmpTree->side-1) * heuristic_evaluation(tmpTree);
                nbContinue++;
                continue;
            }
            compute_attack_squares(tmpTree);
            n_moves = generate_legal_moves(tmpTree, moves);
            if (n_moves==0) {
                tmpResult->score = check(tmpTree) ? -MAX_SCORE : CERTAIN_DRAW;
                nbContinue++;
                continue;
            }
            beg = sizeTree;
            sizeTree += n_moves;
            sum += n_moves;
            if ( (preEvalTrees = realloc(preEvalTrees, sizeTree*sizeof(recTree_t))) == NULL) {
                fprintf(stderr,"realloc error in pre_evaluate()\n");
                exit(1);
            }
            for (int j=beg ; j<sizeTree ; j++) {
                preEvalTrees[j].parentId = i;
                preEvalTrees[j].move = moves[j-beg];
                preEvalTrees[j].tree = malloc(sizeof(tree_t));
                preEvalTrees[j].result = malloc(sizeof(result_t));
                if (!preEvalTrees[j].tree || !preEvalTrees[j].result) {
                    fprintf(stderr,"malloc error in pre_evaluate()\n");
                    exit(1);
                }
                play_move(preEvalTrees[i].tree, moves[j-beg], preEvalTrees[j].tree);
            }
        }
        if (nbContinue == nb_tasks) {
            nb_tasks = 0;
            break;
        }
        nb_tasks = sum;
    }
    beg = sizeTree - nb_tasks;

    /* ------------ TESTING (not parallel) --------------------*/
    for (int i=beg ; i<sizeTree ; i++) {
        evaluate(preEvalTrees[i].tree, preEvalTrees[i].result);
    }
    /*---------------------------------------------------------*/
    // Do parallel thingy
    /* 
     * Master sets up to send tasks to every slave
     * Send & Receive
     * Upon reception, make adjustment (best score thingy)
     */

    /*---------------------------------------------------------*/
    for (int i=sizeTree-1 ;  i>0 ; i--) {
        /* To shorten and clarify expressions */
        recTree_t task = preEvalTrees[i];
        recTree_t parent = preEvalTrees[task.parentId];
        int taskScore = -task.result->score;

        if (taskScore > parent.result->score) {
            parent.result->score = taskScore;
            parent.result->best_move = task.move;
            parent.result->pv_length = task.result->pv_length+1;
            memcpy(parent.result->PV+1, task.result->PV,
                    task.result->pv_length*sizeof(int));
        }
        parent.tree->alpha = MAX(parent.tree->alpha, taskScore);
    }
}

// Only slaves call evaluate
void evaluate(tree_t * T, result_t *result)
{
    node_searched++;

    move_t moves[MAX_MOVES];
    int n_moves;

    result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

    if (test_draw_or_victory(T, result))
        return;

    if (TRANSPOSITION_TABLE && tt_lookup(T, result))     /* la réponse est-elle déjà connue ? */
        return;

    compute_attack_squares(T);

    /* profondeur max atteinte ? si oui, évaluation heuristique */
    if (T->depth == 0) {
        result->score = (2 * T->side - 1) * heuristic_evaluation(T);
        return;
    }

    n_moves = generate_legal_moves(T, &moves[0]);

    /* absence de coups légaux : pat ou mat */
    if (n_moves == 0) {
        result->score = check(T) ? -MAX_SCORE : CERTAIN_DRAW;
        return;
    }

    if (ALPHA_BETA_PRUNING)
        sort_moves(T, n_moves, moves);

    /* évalue récursivement les positions accessibles à partir d'ici */
    for (int i = 0; i < n_moves; i++) {
        tree_t child;
        result_t child_result;

        play_move(T, moves[i], &child);

        evaluate(&child, &child_result);

        int child_score = -child_result.score;

        if (child_score > result->score) { 
            result->score = child_score;
            result->best_move = moves[i];
            result->pv_length = child_result.pv_length + 1;
            for(int j = 0; j < child_result.pv_length; j++)
                result->PV[j+1] = child_result.PV[j];
            result->PV[0] = moves[i];
        }

        if (ALPHA_BETA_PRUNING && child_score >= T->beta)
            break;    

        T->alpha = MAX(T->alpha, child_score);
    }

    if (TRANSPOSITION_TABLE)
        tt_store(T, result);
}


void decide(tree_t * T, result_t *result){
    for (int depth = 1;; depth++) {
        T->depth = depth;
        T->height = 0;
        T->alpha_start = T->alpha = -MAX_SCORE - 1;
        T->beta = MAX_SCORE + 1;

        printf("=====================================\n");
        pre_evaluate(T, result);
        // TODO REMOVE THAT WHEN NOT TESTING
        //if (depth==2) break;

        printf("depth: %d / score: %.2f / best_move : ", T->depth, 0.01 * result->score);
        print_pv(T, result);

        if (DEFINITIVE(result->score))
            break;
    }
}

int main(int argc, char **argv)
{  
    tree_t root;
    result_t result;
    int rank=0, comSize;
    /* MPI Initialization */
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (comSize < 2) {
        fprintf(stderr, "This program requires at least two processes!\n");
        exit(1);
    }


    if (rank==0) {
        if (argc < 2) {
            printf("usage: %s \"4k//4K/4P w\" (or any position in FEN)\n", argv[0]);
            exit(1);
        }

        if (ALPHA_BETA_PRUNING)
            printf("Alpha-beta pruning ENABLED\n");

        if (TRANSPOSITION_TABLE) {
            printf("Transposition table ENABLED\n");
            init_tt();
        }

        parse_FEN(argv[1], &root);
        print_position(&root);


        decide(&root, &result);

        printf("\nDécision de la position: ");
        switch(result.score * (2*root.side - 1)) {
            case MAX_SCORE: printf("blanc gagne\n"); break;
            case CERTAIN_DRAW: printf("partie nulle\n"); break;
            case -MAX_SCORE: printf("noir gagne\n"); break;
            default: printf("BUG\n");
        }

        printf("Node searched: %llu\n", node_searched);

        if (TRANSPOSITION_TABLE)
            free_tt();
    }
    else {
        // TODO
        fprintf(stderr,"Proc %d ne fait rien\n", rank);
    }

    /* MPI Finalization */
    MPI_Finalize();
    return 0;
}
