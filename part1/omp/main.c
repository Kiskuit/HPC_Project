#include "projet.h"
#include <omp.h>

#define OMP_MAX_PROF 1

/* 2017-02-23 : version 1.0 */

unsigned long long int node_searched = 0;

/*void evaluate_omp (tree_t *T, result_t *result, int prof) {
#pragma omp task
{ // OMP BLOCK
if(prof < OMP_PROF_MAX(omp_get_num_threads()))
evaluate_omp (T, result, prof+1);
evaluate (T, result);
} // OMP BLOCK
}*/

void evaluate(tree_t * T, result_t *result, int prof)
{
#pragma omp atomic
    node_searched++;

    move_t moves[MAX_MOVES];
    int n_moves;

    result->score = -MAX_SCORE - 1;
    result->pv_length = 0;

    //printf(">>>>>>>>>>>>>>>>>>>>>>>><<%f\n",result->score);
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
    tree_t child[n_moves];
    result_t child_result[n_moves];
    if(prof < OMP_MAX_PROF){
#pragma omp parallel for
        for (int i = 0; i < n_moves; i++) {
            play_move(T, moves[i], &child[i]);

            evaluate(&child[i], &child_result[i], prof+1);

            int child_score = -child_result[i].score;

#pragma critical CHILD
            { // BLOCK OMP
                if (child_score > result->score) {
                    result->score = child_score;
                    result->best_move = moves[i];
                    result->pv_length = child_result[i].pv_length + 1;
                    for(int j = 0; j < child_result[i].pv_length; j++)
                        result->PV[j+1] = child_result[i].PV[j];
                    result->PV[0] = moves[i];
                }

                // TODO section critique
                T->alpha = MAX(T->alpha, child_score);
            } // BLOCK OMP
        }
    } else {
        for (int i = 0; i < n_moves; i++) {
            play_move(T, moves[i], &child[i]);

            evaluate(&child[i], &child_result[i], prof);

            int child_score = -child_result[i].score;

            if (child_score > result->score) {
                result->score = child_score;
                result->best_move = moves[i];
                result->pv_length = child_result[i].pv_length + 1;
                for(int j = 0; j < child_result[i].pv_length; j++)
                    result->PV[j+1] = child_result[i].PV[j];
                result->PV[0] = moves[i];
            }

            T->alpha = MAX(T->alpha, child_score);
        }

    }
    if (TRANSPOSITION_TABLE)
        tt_store(T, result);
}


void decide(tree_t * T, result_t *result)
{
    for (int depth = 1;; depth++) {
        T->depth = depth;
        T->height = 0;
        T->alpha_start = T->alpha = -MAX_SCORE - 1;
        T->beta = MAX_SCORE + 1;

        printf("=====================================\n");
        //#pragma omp parallel firstprivate(T,result)
        //#pragma omp single
        evaluate(T, result, 0);
        //#pragma omp barrier

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
    return 0;
}
