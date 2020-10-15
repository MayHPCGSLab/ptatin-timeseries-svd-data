
static char help[] = "Load a collection of PETSc vectors and stack them into a tall-skinny matrix\n\n";

#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>

/* Mesh size (129, 33, 257) <xdmf reverse order>, so ni = 257, nk = 129 */
PetscErrorCode SnapshotView(Vec u,const char suffix[])
{
  const PetscInt M = 257;
  const PetscInt N = 33;
  const PetscInt P = 129;
  DM             dm;
  Vec            cu;
  PetscViewer    viewer;
  char           ofile[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,M,N,P,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      1,1,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(dm,0.0,4.0,-1.0,0.0,0.0,2.0);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm,&cu);CHKERRQ(ierr);
  ierr = VecCopy(u,cu);CHKERRQ(ierr);
  
  PetscSNPrintf(ofile,PETSC_MAX_PATH_LEN-1,"%s.vts",suffix);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,ofile,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(cu,viewer);CHKERRQ(ierr);
  
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = VecDestroy(&cu);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SnapshotViewFromFile(const char coor_fname[],const char u_fname[],const char suffix[])
{
  const PetscInt M = 257;
  const PetscInt N = 33;
  const PetscInt P = 129;
  DM             dm;
  Vec            u,coor,_u,_coor;
  PetscViewer    viewer;
  char           ofile[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;
  
  PetscFunctionBegin;
  PetscPrintf(PETSC_COMM_WORLD,"Loading solution: %s\n",u_fname);
  ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,u_fname,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(u,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Loading coordinates: %s\n",coor_fname);
  ierr = VecCreate(PETSC_COMM_WORLD,&coor);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,coor_fname,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(coor,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = DMDACreate3d(PETSC_COMM_WORLD,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,DM_BOUNDARY_NONE,
                      DMDA_STENCIL_BOX,M,N,P,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
                      1,1,NULL,NULL,NULL,&dm);CHKERRQ(ierr);
  ierr = DMSetUp(dm);CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(dm,0.0,1.0,0.0,1.0,0.0,1.0);CHKERRQ(ierr);
  
  ierr = DMCreateGlobalVector(dm,&_u);CHKERRQ(ierr);
  ierr = VecCopy(u,_u);CHKERRQ(ierr);

  ierr = DMGetCoordinates(dm,&_coor);CHKERRQ(ierr);
  ierr = VecCopy(coor,_coor);CHKERRQ(ierr);
  
  PetscSNPrintf(ofile,PETSC_MAX_PATH_LEN-1,"%s.vts",suffix);
  PetscPrintf(PETSC_COMM_WORLD,"Writing output: %s\n",ofile);
  ierr = PetscViewerVTKOpen(PETSC_COMM_WORLD,ofile,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
  ierr = VecView(_u,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = VecDestroy(&_u);CHKERRQ(ierr);
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&coor);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode SnapshotMatCreate(Mat *snapshots)
{
  PetscErrorCode  ierr;
  Vec             u;
  PetscViewer     viewer;
  const PetscInt  steps[] = { 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400 };
  char            fname[PETSC_MAX_PATH_LEN];
  PetscInt        i,len,m,n,k;
  Mat             S;
  PetscInt        *idx;

  PetscFunctionBegin;
  len = sizeof(steps) / sizeof(PetscInt);
  n = len;
  ierr = PetscOptionsGetInt(NULL,NULL,"-truncate",&n,NULL);CHKERRQ(ierr);
  if (n > len) SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_USER,"Truncate value cannot be larger than %D",len);
  
  for (i=0; i<n; i++) {
    PetscSNPrintf(fname,PETSC_MAX_PATH_LEN-1,"data/step%1.6d_energy.pbvec",steps[i]);
    PetscPrintf(PETSC_COMM_WORLD,"[%d] Loading: %s\n",i,fname);
  
    ierr = VecCreate(PETSC_COMM_WORLD,&u);CHKERRQ(ierr);
    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fname,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
    ierr = VecLoad(u,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
    
    if (i == 0) {
      ierr = VecGetSize(u,&m);CHKERRQ(ierr);
      printf("S: m %d, n %d\n",m,n);
      ierr = MatCreate(PETSC_COMM_WORLD,&S);CHKERRQ(ierr);
      ierr = MatSetSizes(S,PETSC_DECIDE,PETSC_DECIDE,m,n);CHKERRQ(ierr);
      ierr = MatSetType(S,MATDENSE);CHKERRQ(ierr);
      ierr = MatSetUp(S);CHKERRQ(ierr);
      
      ierr = PetscCalloc1(m,&idx);CHKERRQ(ierr);
      for (k=0; k<m; k++) { idx[k] = k; }
    }
    {
      const PetscScalar *_u;
      
      ierr = VecGetArrayRead(u,&_u);CHKERRQ(ierr);
      ierr = MatSetValues(S,m,idx,1,&i,_u,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecRestoreArrayRead(u,&_u);CHKERRQ(ierr);
    }
    ierr = MatAssemblyBegin(S,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(S,MAT_FLUSH_ASSEMBLY);CHKERRQ(ierr);
    
    if (i == 0) {
      ierr = SnapshotView(u,"snapshot0");CHKERRQ(ierr);
    }
    
    ierr = VecDestroy(&u);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(S,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(S,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  *snapshots = S;
  
  ierr = PetscFree(idx);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

int main(int argc,char **args)
{
  PetscErrorCode ierr;
  Mat            S;
  PetscInt       step = -1;
  PetscBool      found = PETSC_FALSE;
  
  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  
  ierr = SnapshotMatCreate(&S);CHKERRQ(ierr);

  /* Dump S(i,j) as ascii to stdout */
  /*ierr = MatView(S,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);*/
  
  /* Dump S(i,j) as bindary to file */
  {
    PetscViewer  viewer;
    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"snapshot.pbmat",FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
    ierr = MatView(S,viewer);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);
  }
  
  ierr = PetscOptionsGetInt(NULL,NULL,"-solution_view",&step,&found);CHKERRQ(ierr);
  if (found) {
    char s_fname[PETSC_MAX_PATH_LEN];
    char c_fname[PETSC_MAX_PATH_LEN];
    char ofname[PETSC_MAX_PATH_LEN];
    
    PetscSNPrintf(s_fname,PETSC_MAX_PATH_LEN-1,"data/step%1.6d_energy.pbvec",step);
    PetscSNPrintf(c_fname,PETSC_MAX_PATH_LEN-1,"data/step%1.6d_coor.pbvec",step);
    PetscSNPrintf(ofname,PETSC_MAX_PATH_LEN-1,"step%d_temperature",step);
    PetscPrintf(PETSC_COMM_WORLD,"Snapshot to visualize: %d\n",step);
    ierr = SnapshotViewFromFile(c_fname,s_fname,ofname);CHKERRQ(ierr);
  }
  
  ierr = PetscFinalize();
  return ierr;
}
