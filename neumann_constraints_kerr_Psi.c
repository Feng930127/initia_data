static char help[] = "Einstein constraint equations\n\n";

#include <petsc.h>

#define EPSILON 1.0e-12

typedef struct {
  PetscReal u, v;       //u is conformal factor \psi,
} Field;                //v is stream function \Psi.

typedef struct {
  PetscReal  L,R,       //left and right side of the domain
             U,D,       //upper and down side of the domain
             a,M;       //spin and mass of the black hole
} AppCtx;

/*
We define some functions here, rtuta is the quasi-isotropic radius for
the Kerr spacetime, similarily isotropic radius for the Scharzschild 
spacetime. 
*/
static PetscReal Delta_fun(PetscReal x, PetscReal y,
                           PetscReal a, PetscReal M) {
  PetscReal rtuta, Delta;
  rtuta = x + M + (M*M - a*a)/(4.0*x);
  Delta = a*a - 2.0*M*rtuta + rtuta*rtuta;

  if(PetscAbsReal(Delta) < EPSILON) return(0.);
  else return(Delta);
}

static PetscReal CC_fun(PetscReal x, PetscReal y,
                        PetscReal a, PetscReal M) {
  PetscReal CC;
  CC = (1.0 - (M*M - a*a)/(4*x*x));

  if(PetscAbsReal(CC) < EPSILON) return(0.);
  else return(CC);
}

static PetscReal psi_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, siny, cosy, siny2,
            cosy2, Delta, Sigma, AA, BB, CC;

  rtuta = x + M + (M*M - a*a)/(4.0*x);
  siny  = PetscSinReal(y);
  cosy  = PetscCosReal(y);
  siny2 = PetscPowReal(siny,2.0);
  cosy2 = PetscPowReal(cosy,2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*cosy2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) - a*a*Delta*siny2;
  CC    = CC_fun(x,y,a,M);
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*siny2;

  return PetscPowReal(AA/(x*x*Sigma),1.0/4.0)
                + 0.0*(siny+cosy+cosy2+Delta+Sigma+CC+BB);
}

static PetscReal Psi_exact(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, cosy, siny2,
            cosy2, Sigma;

  rtuta = x + M + (M*M - a*a)/(4.0*x);
  cosy  = PetscCosReal(y);
  siny2 = PetscPowReal(PetscSinReal(y),2.0);
  cosy2 = PetscPowReal(PetscCosReal(y),2.0);
  Sigma = rtuta*rtuta + a*a*cosy2;

  return 2.0*M*a*(cosy2*cosy-3.0*cosy)-2.0*M*a*a*a*cosy*siny2*siny2/Sigma;
}

static PetscReal dNdr(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, siny2, dNdr,
            cosy2, Delta, Sigma, AA, BB, CC;

  rtuta = x + M + (M*M - a*a)/(4.0*x);
  siny2 = PetscPowReal(PetscSinReal(y),2.0);
  cosy2 = PetscPowReal(PetscCosReal(y),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*cosy2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*siny2;
  CC    = CC_fun(x,y,a,M);
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*siny2;

  dNdr = - (AA*AA*CC*(Delta*rtuta + (M - rtuta)*Sigma)
             + a*a*rtuta*(-AA*AA*CC+4.0*AA*CC*M*M*(rtuta*rtuta-Sigma)
             + 2.0*BB*M*M*rtuta*Sigma)*siny2)/(AA*Sigma*(AA*Delta
             - a*a*(AA - 4.0*M*M*rtuta*rtuta)*siny2));

  if(PetscAbsReal(dNdr) < EPSILON) return(0.);
  else return(dNdr);
}

static PetscReal dNdtheta(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, siny, cosy, siny2,
            cosy2, Delta, Sigma, AA, dNdtheta;

  rtuta = x + M + (M*M - a*a)/(4.0*x);
  siny  = PetscSinReal(y);
  cosy  = PetscCosReal(y);
  siny2 = PetscPowReal(PetscSinReal(y),2.0);
  cosy2 = PetscPowReal(PetscCosReal(y),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*cosy2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0) - a*a*Delta*siny2;

  dNdtheta = -a*a*cosy*siny/Sigma + a*a*Delta*cosy*siny/AA;

  if(PetscAbsReal(dNdtheta) < EPSILON) return(0.);
  else return(dNdtheta);
}

static PetscReal qrtheta(PetscReal x, PetscReal y,
                                PetscReal a, PetscReal M){
  PetscReal rtuta, siny, cosy, siny2,
            cosy2, Delta, Sigma, AA, BB, CC;

  rtuta = x + M + (M*M - a*a)/(4.0*x);
  siny  = PetscSinReal(y);
  cosy  = PetscCosReal(y);
  siny2 = PetscPowReal(PetscSinReal(y),2.0);
  cosy2 = PetscPowReal(PetscCosReal(y),2.0);
  Delta = Delta_fun(x,y,a,M);
  Sigma = rtuta*rtuta + a*a*cosy2;
  AA    = PetscPowReal(a*a+rtuta*rtuta,2.0)-a*a*Delta*siny2;
  CC    = (1.0 - (M*M - a*a)/(4*x*x));
  BB    = 4.0*CC*rtuta*(a*a + rtuta*rtuta)
                - a*a*(-2.0*M*CC + 2.0*CC*rtuta)*siny2;

 return BB*(-BB*Sigma*Sigma/(AA*AA) + 4.0*CC*rtuta*Sigma/AA)/(2.0*Sigma*Sigma)
        - (2.0*AA*CC*rtuta*((4.0*CC*rtuta*Sigma)/AA-BB*Sigma*Sigma/(AA*AA)))/(Sigma*Sigma*Sigma)
        + (AA*(-BB*Sigma*Sigma/(AA*AA) + 4.0*CC*rtuta*Sigma/AA))/(2*Sigma*Sigma*x)
        + (1.0/(2.0*Sigma*Sigma))*AA*(-((8.0*BB*CC*Sigma*rtuta)/(AA*AA))
                + (8.0*CC*CC*rtuta*rtuta)/AA + (4.0*CC*CC*Sigma)/AA
                + (2.0*BB*BB*Sigma*Sigma)/(AA*AA*AA)
                + (2.0*(-a*a + M*M)*rtuta*Sigma)/(AA*x*x*x)
                - (Sigma*Sigma*(8.0*CC*CC*rtuta*rtuta + 4.0*CC*CC*(a*a + rtuta*rtuta)
                + (2.0*(-a*a + M*M)*rtuta*(a*a + rtuta*rtuta))/(x*x*x)
                - a*a*(2.0*CC*CC - (M*(-a*a + M*M))/(x*x*x)
                + ((-a*a + M*M)*rtuta)/(x*x*x))*siny2))/(AA*AA))
        + (1.0/(x*x))*((2.0*a*a*AA*cosy*siny*(-((4.0*a*a*Sigma*cosy*siny)/AA)
                + (2.0*a*a*Delta*Sigma*Sigma*cosy*siny)/(AA*AA)))/(Sigma*Sigma*Sigma)
                - (a*a*Delta*cosy*siny*(-((4.0*a*a*Sigma*cosy*siny)/AA)
                + (2.0*a*a*Delta*Sigma*Sigma*cosy*siny)/(AA*AA)))/(Sigma*Sigma)
                + (1.0/(2.0*Sigma*Sigma))*AA*(-((4.0*a*a*Sigma*cosy2)/AA)
                + (2.0*a*a*Delta*Sigma*Sigma*cosy2)/(AA*AA)
                + (4.0*a*a*Sigma*siny2)/AA - (2.0*a*a*Delta*Sigma*Sigma*siny2)/(AA*AA)
                + (8.0*a*a*a*a*cosy2*siny2)/AA - (16.0*a*a*a*a*Delta*Sigma*cosy2*siny2)/(AA*AA)
                + (8.0*a*a*a*a*Delta*Delta*Sigma*Sigma*cosy2*siny2)/(AA*AA*AA)));
}
extern PetscErrorCode InitialState(DM, Vec, AppCtx*);
extern PetscErrorCode ExactSolution(DM, Vec, AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, Field**, Field**, AppCtx*);
extern PetscErrorCode FormJacobianLocal(DMDALocalInfo*, Field**, Mat, Mat, AppCtx*);

//STARTMAIN
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  DM            da;
  SNES          snes;
  AppCtx        user;
  Vec           x, Yexact, err, err1;
  Mat           M1;
  PetscReal     hr, htheta, errnorm, errnorm1;
  DMDALocalInfo info;

  ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;
  user.a  = 0.1;
  user.M  = 1.0;
  user.L  = PetscSqrtReal(user.M*user.M-user.a*user.a)/2;
  user.R  = 8.0;
  user.U  = 0.01;
  user.D  = PETSC_PI-0.01;  
  
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, "Kerr_",
                                "options for Kerr", ""); CHKERRQ(ierr);
  ierr = PetscOptionsReal("-a","The spin of the black hole",
           "jacobian_neumann_constraints_kerr_Psi.c",user.a,&user.a,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-M","The mass of the black hole",
           "jacobian_neumann_constraints_kerr_Psi.c",user.M,&user.M,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-L","The left side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.L,&user.L,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-R","The right side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.R,&user.R,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-U","The upper side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.U,&user.U,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-D","The down side of the domain",
           "jacobian_neumann_constraints_kerr_Psi.c",user.D,&user.D,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd(); CHKERRQ(ierr);

  ierr = DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_STAR,  // for 5-point stencil
               9, 5, PETSC_DECIDE,PETSC_DECIDE,
               2, 1,              // degrees of freedom, stencil width
               NULL, NULL, &da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 0, "u"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 1, "v"); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da, user.L, user.R,
                                user.U, user.D, -1.0, -1.0); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "Solving Einstein constraint equations on %d x %d grid\n", 
			info.mx,info.my); CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&Yexact); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&err); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&err1); CHKERRQ(ierr);
  ierr = InitialState(da,x,&user); CHKERRQ(ierr);
  ierr = ExactSolution(da,Yexact,&user); CHKERRQ(ierr);
  ierr = DMCreateMatrix(da,&M1); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
             (DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  ierr = DMDASNESSetJacobianLocal(da,
             (DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  ierr = SNESComputeJacobianDefault(snes, x, M1, M1, &user); CHKERRQ(ierr);
  //ierr = MatView(M1, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = SNESSolve(snes,NULL,x); CHKERRQ(ierr);
  
  hr = (user.R - user.L) / (PetscReal)(info.mx - 1);
  htheta = (user.D - user.U) / (PetscReal)(info.my - 1);
  ierr = VecWAXPY(err,-1.0,x,Yexact); CHKERRQ(ierr);
  ierr = VecWAXPY(err1,-1.0,x,Yexact); CHKERRQ(ierr);
  ierr = VecScale(err1, hr*htheta); CHKERRQ(ierr);
  ierr = VecNorm(err, NORM_INFINITY, &errnorm); CHKERRQ(ierr);
  ierr = VecNorm(err1, NORM_1, &errnorm1); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "L-Inf(err)=%g, L1(err)=%g\n", errnorm,errnorm1); CHKERRQ(ierr);

  VecDestroy(&x);    VecDestroy(&Yexact); 
  VecDestroy(&err);   VecDestroy(&err1);
  SNESDestroy(&snes);    DMDestroy(&da);
  return PetscFinalize();
}
//ENDMAIN

PetscErrorCode InitialState(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        hr, htheta, x, y, a = user->a, M = user->M;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = VecSet(Y,0.0); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  hr = (user->R - user->L) / (PetscReal)(info.mx - 1);
  htheta = (user->D - user->U) / (PetscReal)(info.my - 1);

  for (j = info.ys; j < info.ys+info.ym; j++) {
      y = j * htheta + user->U;
      for (i = info.xs; i < info.xs+info.xm; i++) {
          x = i * hr + user->L;
          if (i==0||i==info.mx-1||j==0||j==info.my-1) {
                aY[j][i].u = 1.0 + 0.0*psi_exact(x, y, a, M);
                aY[j][i].v = 1.0;//Psi_exact(x, y, a, M);
             }
             else{
                aY[j][i].u = 1.0;//psi_exact(x, y, a, M);
                aY[j][i].v = 1.0;//Psi_exact(x, y, a, M);
             }
      }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}

//Exact Solutions
PetscErrorCode ExactSolution(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        x, y, hr, htheta,
                   a = user->a, M = user->M;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = VecSet(Y,0.0); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  hr = (user->R - user->L) / (PetscReal)(info.mx - 1);
  htheta = (user->D - user->U) / (PetscReal)(info.my - 1);

  for (j = info.ys; j < info.ys+info.ym; j++) {
      y = htheta * j + user->U;
      for (i = info.xs; i < info.xs+info.xm; i++) {
          x = hr * i + user->L;
          aY[j][i].u = psi_exact(x, y, a, M);
          aY[j][i].v = Psi_exact(x, y, a, M);
      }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, Field **aY, 
				Field **aG, AppCtx *user) {
  PetscInt   i, j, mx = info->mx, my = info->my;
  PetscReal  uxx, uyy, u1, u2, u3, vxx, vyy, v1, v2, v3, x, y,
	     L = user->L, R = user->R, U = user->U, D = user->D;
  PetscReal  hr = (user->R - user->L) / (PetscReal)(mx -1),
             htheta = (user->D - user->U) / (PetscReal)(my -1),
             a = user->a, M = user->M, vil, vir, vjl, vjr,
             uil, uir, ujl, ujr;

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = htheta * j + user->U;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hr * i + user->L;
          if (i==0||i==mx-1||j==0||j==my-1) {
              aG[j][i].u = (aY[j][i].u - psi_exact(x, y, a, M));
              aG[j][i].v = (aY[j][i].v - Psi_exact(x, y, a, M));
          }
          else {
              uil = (i==1) ? psi_exact(L, y, a, M) : aY[j][i-1].u;
              uir = (i==mx-2) ? psi_exact(R, y, a, M) : aY[j][i+1].u;
              ujl = (j==1) ? psi_exact(x, U, a, M) : aY[j-1][i].u;
              ujr = (j==my-2) ? psi_exact(x, D, a, M) : aY[j+1][i].u;
              vil = (i==1) ? Psi_exact(L, y, a, M) : aY[j][i-1].v;
              vir = (i==mx-2) ? Psi_exact(R, y, a, M) : aY[j][i+1].v;
              vjl = (j==1) ? Psi_exact(x, U, a, M) : aY[j-1][i].v;
              vjr = (j==my-2) ? Psi_exact(x, D, a, M) : aY[j+1][i].v;

	      uxx = (uil - 2.0 * aY[j][i].u + uir)/(hr)*htheta;
              uyy = (ujl - 2.0 * aY[j][i].u + ujr)/(x*x*htheta)*hr;
              u1 = (uir - uil) / (x)*htheta
                        + (ujr - ujl)/(x*x*2.0*PetscTanReal(y))*hr;
              u2 = 1.0/4.0*qrtheta(x,y,a,M)*aY[j][i].u*hr*htheta;
              u3 = 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (PetscPowReal(1.0/2.0*(vir-vil),2.0)*htheta/hr
		 	+ 1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),2.0)*hr/htheta)
		   *PetscPowReal(aY[j][i].u,-7.0);

              vxx = (vil - 2.0 * aY[j][i].v + vir)/(hr)*htheta;
              vyy = (vjl - 2.0 * aY[j][i].v + vjr)/(x*x*htheta)*hr;
              v1 = - 2.0/x*1.0/(2.0)*(vir - vil)*htheta
                        - 3.0/(x*x)*1.0/PetscTanReal(y)*1.0/(2.0)*(vjr - vjl)*hr;
              v2 = dNdr(x,y,a,M)*1.0/(2.0)*(vir - vil)*htheta
                        + dNdtheta(x,y,a,M)*1.0/(2.0*x*x)*(vjr - vjl)*hr;
              v3 = - 6.0/(aY[j][i].u)*1.0/(2.0)*(uir - uil)
                                *1.0/(2.0*hr)*(vir - vil)*htheta
                   - 6.0/(aY[j][i].u)*1.0/(2.0*x*x)*(ujr - ujl)
                                *1.0/(2.0*htheta)*(vjr - vjl)*hr;
	      aG[j][i].u = ((uxx + uyy) + u1 + u2 + u3);
              aG[j][i].v = ((vxx + vyy) + v1 + v2 + v3);
          }
      }
    }
    return 0;
}

PetscErrorCode FormJacobianLocal(DMDALocalInfo *info, Field **aY,
                                    Mat J, Mat P, AppCtx *user) {
  PetscErrorCode ierr;
  PetscInt    i, j, mx = info->mx, my = info->my;
  PetscReal   v[5],L=user->L, R=user->R, U=user->U, D=user->D;
  MatStencil  col[5],row;
  PetscReal   hr = (user->R - user->L) / (PetscReal)(mx -1),
              htheta = (user->D - user->U) / (PetscReal)(my -1),
              a = user->a, M = user->M, vil, vir, vjl, vjr,
              uil, uir, ujl, ujr, x, y;
  
  for (j = info->ys; j < info->ys+info->ym; j++) {
      row.j = j;  y = htheta * j + user->U;
      for (i = info->xs; i < info->xs+info->xm; i++) {
          row.i = i;  x = hr * i + user->L;
	  if (i==0||i==mx-1||j==0||j==my-1) {
	      row.c = 0;  col[0].c = 0; 
              col[0].i = i; col[0].j = j; v[0] = 1.0;
	      ierr = MatSetValuesStencil(P,1,&row,1,col,v,INSERT_VALUES); CHKERRQ(ierr);
	      row.c = 1;  col[0].c = 1; 
              col[0].i = i; col[0].j = j; v[0] = 1.0;
              ierr = MatSetValuesStencil(P,1,&row,1,col,v,INSERT_VALUES); CHKERRQ(ierr);
          }
	  else {
  	  uil = (i==1) ? psi_exact(L, y, a, M) : aY[j][i-1].u;
          uir = (i==mx-2) ? psi_exact(R, y, a, M) : aY[j][i+1].u;
          ujl = (j==1) ? psi_exact(x, U, a, M) : aY[j-1][i].u;
          ujr = (j==my-2) ? psi_exact(x, D, a, M) : aY[j+1][i].u;
          vil = (i==1) ? Psi_exact(L, y, a, M) : aY[j][i-1].v;
          vir = (i==mx-2) ? Psi_exact(R, y, a, M) : aY[j][i+1].v;
          vjl = (j==1) ? Psi_exact(x, U, a, M) : aY[j-1][i].v;
          vjr = (j==my-2) ? Psi_exact(x, D, a, M) : aY[j+1][i].v;
	  //psi equation
	  row.c = 0;  col[0].c = 0;  col[1].c = 0;
	  col[2].c = 0;	col[3].c = 0; col[4].c = 0;
          col[0].i = i; col[0].j = j; 
	  v[0] = - 2.0*(htheta/hr + hr/htheta/(x*x)) + 1.0/4.0*qrtheta(x,y,a,M)*hr*htheta
	  	    + 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                       * (PetscPowReal(1.0/(2.0)*(vir-vil),2.0)*htheta/hr
                           + 1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),2.0)*hr/htheta)
                       * PetscPowReal(aY[j][i].u,-8.0)*(-7.0);
	  col[1].i = i-1; col[1].j = j; 
	  v[1] = (i>1) ? 1.0*htheta/hr - 1.0/(x)*htheta : 0.0; 
	  col[2].i = i+1; col[2].j = j; 
	  v[2] = (i<mx-2) ? 1.0*htheta/hr + 1.0/(x)*htheta : 0.0;
   	  col[3].i = i; col[3].j = j-1; 
	  v[3] = (j>1) ? 1.0/(x*x)*hr/htheta - 1.0/(2.0*x*x)*hr/PetscTanReal(y) : 0.0;
   	  col[4].i = i; col[4].j = j+1; 
	  v[4] = (j<my-2) ? 1.0/(x*x)*hr/htheta + 1.0/(2.0*x*x)*hr/PetscTanReal(y) : 0.0;
   	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
	  
	  row.c = 0;  col[0].c = 1;  col[1].c = 1;
          col[2].c = 1; col[3].c = 1; col[4].c = 1;
	  col[0].i = i; col[0].j = j;
 	  v[0]= 0.0;
	  col[1].i = i-1; col[1].j = j;
	  v[1] = (i>1) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (-1.0*PetscPowReal(1.0/2.0*(vir-vil),1.0)*htheta/hr)
                   	* PetscPowReal(aY[j][i].u,-7.0) : 0.0;
	  col[2].i = i+1; col[2].j = j;
          v[2] = (i<mx-2) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (1.0*PetscPowReal(1.0/2.0*(vir-vil),1.0)*htheta/hr)
                        * PetscPowReal(aY[j][i].u,-7.0) : 0.0;
 	  col[3].i = i; col[3].j = j-1;
   	  v[3] = (j>1) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (- 1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),1.0)*hr/htheta)
                   *PetscPowReal(aY[j][i].u,-7.0) : 0.0;
	  col[4].i = i; col[4].j = j+1;
          v[4] = (j<my-2) ? 1.0/16.0*PetscPowReal(x*PetscSinReal(y),-4.0)
                   * (1.0/(x*x)*PetscPowReal(1.0/(2.0)*(vjr-vjl),1.0)*hr/htheta)
                   *PetscPowReal(aY[j][i].u,-7.0) : 0.0;
	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);


 	  //Y equation
	  row.c = 1;  col[0].c = 1;  col[1].c = 1;
          col[2].c = 1; col[3].c = 1; col[4].c = 1;
          col[0].i = i; col[0].j = j; 
	  v[0] = - 2.0*(htheta/hr + hr/htheta/(x*x));
 	  col[1].i = i-1; col[1].j = j; 
	  v[1] = (i>1) ? 1.0*htheta/hr + 2.0/(2.0*x)*htheta - dNdr(x,y,a,M)/2.0*htheta 
			+ 6.0/(2.0*aY[j][i].u)*htheta/(2.0*hr)*(uir-uil) : 0.0;
	  col[2].i = i+1; col[2].j = j;
	  v[2] = (i<mx-2) ? 1.0*htheta/hr - 2.0/(2.0*x)*htheta + dNdr(x,y,a,M)/2.0*htheta
			- 6.0/(2.0*aY[j][i].u)*htheta/(2.0*hr)*(uir-uil) : 0.0;
	  col[3].i = i; col[3].j = j-1;
  	  v[3] = (j>1) ? 1.0/(x*x)*hr/htheta + 3.0/(2.0*x*x*PetscTanReal(y))*hr
			-dNdtheta(x,y,a,M)/(2.0*x*x)*hr
			+ 6.0/(2.0*aY[j][i].u*x*x)*hr/(2.0*htheta)*(ujr-ujl) : 0.0;
	  col[4].i = i; col[4].j = j+1;
	  v[4] = (j<my-2) ? 1.0/(x*x)*hr/htheta - 3.0/(2.0*x*x*PetscTanReal(y))*hr
                	+ dNdtheta(x,y,a,M)/(2.0*x*x)*hr 
			- 6.0/(2.0*aY[j][i].u*x*x)*hr/(2.0*htheta)*(ujr-ujl) : 0.0;
	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);

	  row.c = 1;  col[0].c = 0;  col[1].c = 0;
          col[2].c = 0; col[3].c = 0; col[4].c = 0;
          col[0].i = i; col[0].j = j;
	  v[0] = 6.0/(aY[j][i].u*aY[j][i].u)*1.0/(2.0)*(uir-uil)
                        *1.0/(2.0*hr)*(vir-vil)*htheta
                + 6.0/(aY[j][i].u*aY[j][i].u)*1.0/(2.0*x*x)*(ujr-ujl)
                        *1.0/(2.0*htheta)*(vjr-vjl)*hr;
	  col[1].i = i-1; col[1].j = j;
	  v[1] = (i>1) ? - 6.0/(aY[j][i].u)*1.0/(2.0)*(-1.0)
                                *1.0/(2.0*hr)*(vir - vil)*htheta : 0.0;
	  col[2].i = i+1; col[2].j = j;
	  v[2] = (i<mx-2) ? - 6.0/(aY[j][i].u)*1.0/(2.0)*(1.0)
                                *1.0/(2.0*hr)*(vir - vil)*htheta : 0.0;
	  col[3].i = i; col[3].j = j-1;
	  v[3] = (j>1) ? - 6.0/(aY[j][i].u)*1.0/(2.0*x*x)*(-1.0)
                                *1.0/(2.0*htheta)*(vjr - vjl)*hr : 0.0;
	  col[3].i = i; col[4].j = j+1;
 	  v[4] = (j<my-2) ? - 6.0/(aY[j][i].u)*1.0/(2.0*x*x)*(1.0)
                                *1.0/(2.0*htheta)*(vjr - vjl)*hr : 0.0;
	  ierr = MatSetValuesStencil(P,1,&row,5,col,v,INSERT_VALUES); CHKERRQ(ierr);
	  }
        }
    }
  ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  //ierr = MatView(P,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  if (J != P) {
      ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
      ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  }
  return 0;
}
