
static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol : use a random exact solution vector\n\
  -view_exact_sol   : write exact solution vector to stdout\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";

#include <petscksp.h>

#include <mpi.h>



#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x,b,w;  /* approx solution, RHS, exact solution */
  Mat            A;        /* linear system matrix */
  KSP            ksp;     /* linear solver context */
  PetscInt       i,Istart,Iend,its,N;
  PetscErrorCode ierr;
#if defined(PETSC_USE_LOG)
  PetscLogStage stage;
#endif

  PetscInitialize(&argc,&args,(char*)0,help);
//  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
//  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.

     Performance tuning note:  For problems of substantial size,
     preallocation of matrix memory is crucial for attaining good
     performance. See the matrix chapter of the users manual for details.
  */

int myid,size;
MPI_Comm_rank(MPI_COMM_WORLD, &myid);
MPI_Comm_size(MPI_COMM_WORLD, &size);
int aa = 1;
while(aa==0)
{

}
FILE *fp;
fp = fopen("./data/matrix.dat","r");
fscanf(fp,"%d\n",&N);
int myi,myj;
double val;
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
ierr = MatSetFromOptions(A);CHKERRQ(ierr);
if(myid==0){
	ierr = MatSetSizes(A,N/size+N%size,N/size+N%size,N,N);CHKERRQ(ierr);
}
else{
	ierr = MatSetSizes(A,N/size,N/size,N,N);CHKERRQ(ierr);
}

//	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,N,N);CHKERRQ(ierr);



  /*
     Currently, all PETSc parallel matrix formats are partitioned by
     contiguous chunks of rows across the processors.  Determine which
     rows of the matrix are locally owned.
  */
ierr = MatMPIAIJSetPreallocation(A,8,NULL,7,NULL); CHKERRQ(ierr); 
ierr = MatSeqAIJSetPreallocation(A,8,NULL); CHKERRQ(ierr); 




  ierr = MatGetOwnershipRange(A,&Istart,&Iend);CHKERRQ(ierr);

  /*
     Set matrix elements for the 2-D, five-point stencil in parallel.
      - Each processor needs to insert only elements that it owns
        locally (but any non-local elements will be sent to the
        appropriate processor during matrix assembly).
      - Always specify global rows and columns of matrix entries.

     Note: this uses the less common natural ordering that orders first
     all the unknowns for x = h then for x = 2h etc; Hence you see J = Ii +- n
     instead of J = I +- m as you might expect. The more standard ordering
     would first do all variables for y = h, then y = 2h etc.

   */
  ierr = PetscLogStageRegister("Assembly", &stage);CHKERRQ(ierr);
  ierr = PetscLogStagePush(stage);CHKERRQ(ierr);

while(fscanf(fp,"(%d,%d) %le\n",&myi,&myj,&val)!=EOF){
	if(myi<Iend && myi>= Istart){
		ierr = MatSetValues(A,1,&myi,1,&myj,&val,ADD_VALUES);CHKERRQ(ierr);
	}
}
 
fclose(fp); 
 
 
 
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = PetscLogStagePop();CHKERRQ(ierr);


  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  
  if(myid==0){
  ierr = VecSetSizes(x,N/size+N%size,N);CHKERRQ(ierr);
  }
  else{
  ierr = VecSetSizes(x,N/size,N);CHKERRQ(ierr);
  }

//  ierr = VecSetSizes(x,PETSC_DECIDE,N);CHKERRQ(ierr);



ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);
	
ierr = VecGetOwnershipRange(b,&Istart,&Iend);CHKERRQ(ierr);;
fp = fopen("./data/rhs.dat","r");
myi = 0;

fscanf(fp,"%d",&N);
while(fscanf(fp,"%le\n",&val)!=EOF){
	if(myi<Iend && myi>= Istart){
		ierr = VecSetValues(b,1,&myi,&val,ADD_VALUES);CHKERRQ(ierr);
	}
	++myi;
}
ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
fclose(fp);





  /*
     Create linear solver context
  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */
  ierr = KSPSetTolerances(ksp,1.e-20,1.e-40,PETSC_DEFAULT,
                          PETSC_DEFAULT);CHKERRQ(ierr);

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Check the error
  */
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  /*
     Print convergence information.  PetscPrintf() produces a single
     print statement from all processes that share a communicator.
     An alternative is PetscFPrintf(), which prints to a file.
  */




VecScatter ctx;

ierr = VecCreateSeq(PETSC_COMM_SELF,N,&w);CHKERRQ(ierr);
VecScatterCreateToZero(x,&ctx,&w);
VecScatterBegin(ctx,x,w,INSERT_VALUES,SCATTER_FORWARD);	
VecScatterEnd(ctx,x,w,INSERT_VALUES,SCATTER_FORWARD);	
VecScatterDestroy(&ctx);

printf("scatteredi\n");
//VecView(w, PETSC_VIEWER_STDOUT_WORLD);

//MatView(A, PETSC_VIEWER_STDOUT_WORLD);

PetscScalar *array;
ierr = VecGetArray(w,&array);CHKERRQ(ierr);;

FILE* outfp,*fp1;
double myx,myy;



if(myid==0){
outfp = fopen("out.dat","w");
fp1 = fopen("./data/point.dat","r");
fscanf(fp1,"%d\n",&i);
for(i=0;i<N;++i){
	fscanf(fp1,"%le %le\n",&myx,&myy);
	if(i==0)
		printf("%le %le %le\n",myx,myy,array[i]);
	ierr = PetscFPrintf(PETSC_COMM_WORLD,outfp,"%le %le %le\n",myx,myy,array[i]);CHKERRQ(ierr);
}

fclose(fp1);
fclose(outfp);
}
ierr = VecRestoreArray(w,&array);CHKERRQ(ierr);
  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
ierr = VecDestroy(&x);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);  ierr = MatDestroy(&A);CHKERRQ(ierr);
 
ierr = VecDestroy(&w);CHKERRQ(ierr);
  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary).
  */
  ierr = PetscFinalize();
  return 0;
}
