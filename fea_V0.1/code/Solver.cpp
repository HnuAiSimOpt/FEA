
//
#include"Solver.h"

void SuperLUsolver(SpMatrix& K, DoubleArray& F, DoubleArray& dis)
{
	if (F.SquareNorm() == 0)
	{
		printf("��������\n");
	}
	else
	{
		SuperMatrix A, L, U, B;
		int* perm_r;     /* row permutations from partial pivoting */
		int* perm_c;     /* column permutation vector */
		superlu_options_t options;
		SuperLUStat_t stat;

		/* Initialize matrix A. */
		int* ia = K.GetColPtr().GetPointer();
		int* ja = K.GetRowIndex().GetPointer();
		double* a = K.GetNzValues().GetPointer();
		double* rhs = F.GetPointer();
		int neqs = F.getSize();
		int nz = K.GetNZ();

		/*Create matrix A in the format expected by SuperLU.*/
		dCreate_CompCol_Matrix(&A, neqs, neqs, nz, a, ja, ia, SLU_NC, SLU_D, SLU_GE);

		/* Create right-hand side matrix B. */
		int nrhs = 1;
		dCreate_Dense_Matrix(&B, neqs, nrhs, rhs, neqs, SLU_DN, SLU_D, SLU_GE);
		if (!(perm_r = intMalloc(neqs))) ABORT("Malloc fails for perm_r[].");
		if (!(perm_c = intMalloc(neqs))) ABORT("Malloc fails for perm_c[].");

		/* Set the default input options. */
		set_default_options(&options);
		options.ColPerm = COLAMD;
		//  ����Ԫ�������ʱ��options.ColPerm����ΪCOLAMD�ٶȿ�����ߺü���(SuperLUĬ�ϵ�ColPermֵ����COLAMD)��
		//	ColPer������ug�е������ǣ�Specifies how to permute the columns of the matrix for sparsity preservation.
		//	��Ҫ�漰A����ϡ�财��ı任��ʽ�����������Ҫ����A�����������ѡ���漰����֪ʶ��С�ܲ��ţ�ʵ�ڲ�����

		/* Initialize the statistics variables. */
		StatInit(&stat);

		/* Solve the linear system. */
		int info;
		dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

		dis = F;
		/* De-allocate storage */
		//SUPERLU_FREE(rhs);
		SUPERLU_FREE(perm_r);
		SUPERLU_FREE(perm_c);
		//Destroy_CompCol_Matrix(&A);
		//Destroy_SuperMatrix_Store(&B);
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		StatFree(&stat);
	}
}