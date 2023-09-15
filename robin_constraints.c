static char help[] = "Einstein constraint equations\n\n";

#include <petsc.h>

#define EPSILON 1.0e-12

typedef struct {
  PetscReal u, v;       //u is conformal factor \psi,
} Field;                //v is stream function \Psi.

typedef struct {
  PetscReal  L,R,       //left and right side of the domain
             U,D,       //upper and down side of the domain
             a,M,
             sig1,sig2,
             coer,coethe;       //spin and mass of the black hole
} AppCtx;

static PetscReal barrier(PetscReal x, PetscReal y, 
                                PetscReal rhori) {
  PetscReal barrier, a, r1,r2;
  
  a = PetscPowReal(10.0,20.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y - 3)*(y - 3)) - rhori;
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y + 3)*(y + 3)) - rhori;
  barrier = - 1.5/PETSC_PI*PetscAtanReal(r1*a)
            - 1.5/PETSC_PI*PetscAtanReal(r2*a) + 1.5/2 + 1.5/2;
  if(PetscAbsReal(barrier) < 0.1) return(0.);
  else return(barrier);
}

static PetscReal barrier_rup(PetscReal x, PetscReal y, 
			   PetscReal r0, PetscReal epsilon) {
  PetscReal barrier_rup, a, r1, r2;

  a = PetscPowReal(10.0, 20.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y - 3)*(y - 3)) 
		- (r0 - epsilon);
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y - 3)*(y - 3)) 
		- (r0 + epsilon);
  barrier_rup = 1.0/PETSC_PI*(PetscAtanReal(a*r1)
		- PetscAtanReal(a*r2)); 
  
  return(barrier_rup);
}

static PetscReal barrier_rdown(PetscReal x, PetscReal y,
                             PetscReal r0, PetscReal epsilon) {
  PetscReal barrier_rdown, a, r1, r2;

  a = PetscPowReal(10.0, 20.0);
  r1 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y + 3)*(y + 3)) 
		- (r0 - epsilon);
  r2 = PetscSqrtReal((x - 0.001)*(x - 0.001) + (y + 3)*(y + 3)) 
		- (r0 + epsilon);
  barrier_rdown = 1.0/PETSC_PI*(PetscAtanReal(a*r1)
                  - PetscAtanReal(a*r2));
  
  return(barrier_rdown);
} 

static PetscReal barrier_thetaup(PetscReal x, PetscReal y,
			       PetscReal theta0, PetscReal epsilon) {
  PetscReal barrier_thetaup, a, theta1, theta2;
  
  a = PetscPowReal(10.0, 20.0);
  theta1 = PetscAtanReal((y - 3.0)/x) - (theta0 - epsilon);
  theta2 = PetscAtanReal((y - 3.0)/x) - (theta0 + epsilon);
  barrier_thetaup = 1.0/PETSC_PI*(PetscAtanReal(a*theta1)
                - PetscAtanReal(a*theta2));
  
  return (barrier_thetaup);
}

static PetscReal barrier_thetadown(PetscReal x, PetscReal y,
                               PetscReal theta0, PetscReal epsilon) {
  PetscReal barrier_thetadown, a, theta1, theta2;

  a = PetscPowReal(10.0, 20.0);
  theta1 = PetscAtanReal((y + 3.0)/x) - (theta0 - epsilon);
  theta2 = PetscAtanReal((y + 3.0)/x) - (theta0 + epsilon);
  barrier_thetadown = 1.0/PETSC_PI*(PetscAtanReal(a*theta1)
                - PetscAtanReal(a*theta2));
  
  return (barrier_thetadown);
} 

static PetscReal qro2z2(PetscReal ro, PetscReal z,
                        PetscReal a, PetscReal sig1,
                        PetscReal sig2) {
  PetscReal exp1, part1, part2, part3, part4, qroz;

  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  part1 = -6.0*a*ro*exp1/(sig1*sig1);
  part2 = -2.0*a*ro*exp1/(sig2*sig2);
  part3 = 4.0*a*ro*z*z*exp1/(sig2*sig2*sig2*sig2);
  part4 = 4.0*a*ro*ro*ro*exp1/(sig1*sig1*sig1*sig1);

  qroz = part1 + part2 + part3 + part4;
  if(PetscAbsReal(qroz) < EPSILON) return(0.);
  else return(qroz);
}

static PetscReal qro(PetscReal ro, PetscReal z,
		     PetscReal a, PetscReal sig1,
		     PetscReal sig2) {
  PetscReal qro, exp1;
  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  qro = a*exp1 - 2.0*a*ro*ro*exp1/(sig1*sig1);

  return (qro);
}

static PetscReal qz(PetscReal ro, PetscReal z,
		    PetscReal a, PetscReal sig1,
		    PetscReal sig2) {
  PetscReal qz, exp1;
  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  qz = -2.0*a*ro*z*exp1/(sig2*sig2);

  return (qz);
}

static PetscReal q_func(PetscReal ro, PetscReal z,
		 	PetscReal a, PetscReal sig1,
			PetscReal sig2) {
  PetscReal exp1, qq;
  exp1 = PetscExpReal(-ro*ro/(sig1*sig1) - z*z/(sig2*sig2));
  qq = a*ro*exp1;
  
  return (qq);
}
 
static PetscReal NN(PetscReal ro, PetscReal z) {
  PetscReal NN;

  NN = 1.0;
  return NN;
}

static PetscReal dNdro(PetscReal ro, PetscReal z) {
  PetscReal dNdro;

  dNdro = 0.0;
  return dNdro;
}

static PetscReal dNdz(PetscReal ro, PetscReal z) {
  PetscReal dNdz;

  dNdz = 0.0;
  return dNdz;
}
 
extern PetscErrorCode InitialState(DM, Vec,  AppCtx*);
extern PetscErrorCode FormFunctionLocal(DMDALocalInfo*, Field**,
                                           Field**,  AppCtx*);
extern PetscErrorCode indexfunc_up(DMDALocalInfo*, PetscReal, PetscReal,PetscReal,
				PetscReal*, Field**,  AppCtx*); 
extern PetscErrorCode indexfunc_down(DMDALocalInfo*, PetscReal, PetscReal,PetscReal,
                                PetscReal*, Field**,  AppCtx*);
extern PetscErrorCode AngMom_ADMMass(DMDALocalInfo*, DM, Vec, PetscReal, AppCtx*);


//STARTMAIN
int main(int argc,char **argv) {
  PetscErrorCode ierr;
  DM            da;
  SNES          snes;
  AppCtx        user;
  Vec           x, Yexact, err, err1, y;
  PetscReal     AngMom = 0.0, hx,hy, errnorm, errnorm1, errnorm2;
  DMDALocalInfo info;
  Field   **aY;
  PetscInt i,j;

  ierr = PetscInitialize(&argc,&argv,NULL,help); if (ierr) return ierr;
  user.a  = 4.0;
  user.M  = 1.0;
  user.L  = 0.001;
  user.R  = 20.0;
  user.U  = 10.0;
  user.D  = -10.0;
  user.sig1 = 1.0;
  user.sig2 = 1.0;
  user.coer = 0.7;
  user.coethe = 0.5;
 
  ierr = DMDACreate2d(PETSC_COMM_WORLD,
               DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
               DMDA_STENCIL_STAR,  // for 5-point stencil
               11, 11, PETSC_DECIDE,PETSC_DECIDE,
               2, 1,              // degrees of freedom, stencil width
               NULL, NULL, &da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 0, "psi"); CHKERRQ(ierr);
  ierr = DMDASetFieldName(da, 1, "Omega"); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DMDASetUniformCoordinates(da, user.L, user.R,
                                user.D, user.U, -1.0, -1.0); CHKERRQ(ierr);
  ierr = DMSetApplicationContext(da,&user); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "Solving Einstein constraint equations on %d x %d grid\n", 
			info.mx,info.my); CHKERRQ(ierr);

  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&x); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&Yexact); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&err); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&err1); CHKERRQ(ierr);
  ierr = VecDuplicate(x,&y); CHKERRQ(ierr);
  ierr = InitialState(da,x,&user); CHKERRQ(ierr);

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes); CHKERRQ(ierr);
  ierr = SNESSetDM(snes,da); CHKERRQ(ierr);
  ierr = DMDASNESSetFunctionLocal(da,INSERT_VALUES,
             (DMDASNESFunction)FormFunctionLocal,&user); CHKERRQ(ierr);
  //ierr = DMDASNESSetJacobianLocal(da,
    //         (DMDASNESJacobian)FormJacobianLocal,&user); CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
  
  ierr = SNESSolve(snes,NULL,x); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,x,&aY); CHKERRQ(ierr);

  ierr = SNESComputeFunction(snes, x, y); CHKERRQ(ierr);
  ierr = AngMom_ADMMass(&info, da, x, AngMom, &user);
 
  hx = (user.R - user.L) / (PetscReal)(info.mx - 1);
  hy = (user.U - user.D) / (PetscReal)(info.my - 1);
  ierr = VecNorm(y, NORM_INFINITY, &errnorm); CHKERRQ(ierr);
  ierr = VecNorm(y, NORM_1, &errnorm1); CHKERRQ(ierr);

  ierr = VecNorm(y, NORM_INFINITY, &errnorm2); CHKERRQ(ierr);
for (j = info.ys; j < info.ys+info.ym; j++) {
       for (i = info.xs; i < info.xs + 1; i++) {
           //PetscPrintf(PETSC_COMM_WORLD, "psi[%d][%d] = %10f,   Omega[%d][%d] = %10f\n", j,i,aY[j][i].u,j,i,aY[j][i].v);
       }
   }
  ierr = PetscPrintf(PETSC_COMM_WORLD,
           "L-Inf(res)=%g, L1(res)=%g, L-Inf(res)=%g,Angmom = \n", 
                    errnorm,errnorm1*hx*hy,errnorm2, AngMom); CHKERRQ(ierr);

  VecDestroy(&x);    VecDestroy(&Yexact); 
  VecDestroy(&err);   VecDestroy(&err1);
  SNESDestroy(&snes);    DMDestroy(&da);
  return PetscFinalize();
}
//ENDMAIN

//Initial state
PetscErrorCode InitialState(DM da, Vec Y, AppCtx* user) {
  PetscErrorCode   ierr;
  DMDALocalInfo    info;
  PetscInt         i,j;
  PetscReal        x, y, hx, hy, r0, rup, rdown;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = VecSet(Y, 0.0); CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da, &info); CHKERRQ(ierr);
  ierr = DMDAGetCoordinateArray(da, &aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da, Y, &aY); CHKERRQ(ierr);
  hx = (user->R - user->L) / (PetscReal)(info.mx - 1);
  hy = (user->U - user->D) / (PetscReal)(info.my - 1);

  for (j = info.ys; j < info.ys+info.ym; j++) {
      y = j * hy + user->D;
      for (i = info.xs; i < info.xs+info.xm; i++) {
          x = i * hx + user->L;
          if (i == info.mx - 1 || j == 0 || j == info.my - 1) {
              r0 = PetscSqrtReal(x*x + y*y);
              aY[j][i].u = 1.0;//x*x + y*y;
		      aY[j][i].v = 1.0/(r0*r0*r0);
          }
          else {
              rup = PetscSqrtReal(((x - 0.001)*(x - 0.001)) + (y - 3)*(y - 3));
              rdown = PetscSqrtReal(((x - 0.001)*(x - 0.001)) + (y + 3)*(y + 3));
              if (barrier(x, y, 1.5 - 0.52*hx) < 0.2){
                   aY[j][i].u = (1.0 + 1.2/(2.0*rup) + 1.2/(2.0*rdown))*(1.5 - barrier(x, y, 1.5 - 0.52*hx))/1.5
                               + 1.8*barrier(x, y, 1.5 - 0.52*hx);
		           aY[j][i].v =  0.02*barrier(x, y, 1.5 - 0.52*hx);

              }
              else {
                  aY[j][i].u = 1.8;
                  aY[j][i].v = 0.02*barrier(x, y, 1.5 - 0.52*hx);
              }
          }
      }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}

//index function for y > 0 part.
PetscErrorCode indexfunc_up(DMDALocalInfo *info, PetscReal theta,
                        PetscReal r1, PetscReal r2, PetscReal *value, 
  			Field **aY, AppCtx *user) {
  PetscInt   i, j, mx = info->mx, my = info->my, k = 0, l = 0;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx - 1),
             hy = (user->U - user->D) / (PetscReal)(my - 1),
             x, y, coer = user->coer, coethe = user->coethe, 
             barrier1, barrier2;

  for (j = info->ys; j < info->ys + info->ym; j++) {
    y = hy * j + user->D;
    for (i = info->xs; i < info->xs + info->xm; i++) {
        x = hx * i + user->L; 
        barrier1 = barrier_rup(x, y, r1, coer*hx)
                    *barrier_thetaup(x, y, theta, coethe*hx);
 	    barrier2 = barrier_rup(x, y, r2, coer*hx)
                    *barrier_thetaup(x, y, theta, coethe*hx);
	    if (barrier1 > 0.9) {
            value[0] = aY[j][i].u;
            value[2] = aY[j][i].v;
		    k = k + 1;
        }
	    else if (barrier2 > 0.9) {
		    value[1] = aY[j][i].u;
            value[3] = aY[j][i].v;
		    l = l + 1;
	    }
	    else if (k == 1 && l == 1) {
	        goto here;
	    }
    }
  }
  here:
  if (l == 0 || k == 0) {
      PetscPrintf(PETSC_COMM_WORLD, "l = %g,k = %g\n", l,k);
  }
  return 0.;
}

//index function for y < 0 part.
PetscErrorCode indexfunc_down(DMDALocalInfo *info, PetscReal theta,
                        PetscReal r1, PetscReal r2, PetscReal *value, 
			Field **aY, AppCtx *user) {
  PetscInt   i, j, mx = info->mx, my = info->my, k = 0, l = 0;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx - 1),
             hy = (user->U - user->D) / (PetscReal)(my - 1),
             x, y, coer = user->coer, coethe = user->coethe, 
             barrier1, barrier2;

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hx * i + user->L; 
          barrier1 = barrier_rdown(x, y, r1, coer*hx)
                        *barrier_thetadown(x, y, theta, coethe*hx);
	      barrier2 = barrier_rdown(x, y, r2, coer*hx)
                        *barrier_thetadown(x, y, theta, coethe*hx);
          if (barrier1 > 0.9) {
                value[0] = aY[j][i].u;
                value[2] = aY[j][i].v;
		        k = k + 1;
          }     
	      else if (barrier2 > 0.9) {
		        value[1] = aY[j][i].u;
                value[3] = aY[j][i].v;
		        l = l + 1;
	      }
	      else if (k == 1 && l == 1) {
		        goto here;
	      }
       }   
  }   
  here:
  if (l == 0 || k == 0){
    PetscPrintf(PETSC_COMM_WORLD, "l = %g,k = %g\n", l,k);
  }
  return 0.;
}

PetscErrorCode FormFunctionLocal(DMDALocalInfo *info, Field **aY, Field **aG, AppCtx *user) {
  PetscErrorCode   ierr;
  PetscInt   i, j, mx = info->mx, my = info->my, n, k;
  PetscReal  x, y, barrier_rup1, value[4], uxx, uyy,
 	         uxl, vxl, vxx,vyy, barrier_rdown1, barrier_thetaup1,
	         barrier_thetadown1, u1, u2, u3, v1, v2, v3,
	         uxr, uyl, uyr, vxr, vyl, vyr, part, theta = 0.0,
	         a = user->a, sig1 = user->sig1, sig2 = user->sig2,
             coer = user->coer, dOmega, AngMom = 0.0;
  PetscReal  hx = (user->R - user->L) / (PetscReal)(mx - 1),
	         hy = (user->U - user->D) / (PetscReal)(my - 1);

  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
	      x = hx * i + user->L;
	      barrier_rup1 = barrier_rup(x, y, 1.5 - 2.0*coer*hx, coer*hx);
          barrier_rdown1 = barrier_rdown(x, y, 1.5 - 2.0*coer*hx, coer*hx);
	      if (barrier_rup1 > 0.9) { 
              theta = PetscAtanReal((y - 3.0)/x);
		      ierr = indexfunc_up(info, theta, 1.5 + 2.0*coer*hx,
				        1.5, value, aY, user); CHKERRQ(ierr);
		      part = 1.0/(value[1]*PetscExpReal(q_func(x,y,a,sig1,sig2))*1.5);
		      //aY[j][i].u = aY[j][i].u;
              aY[j][i].u = value[0] + 4.0/3.0*coer*hx*part*(x*qro(x,y,a,sig1,sig2)
			  	        + (y - 3.0)*qz(x,y,a,sig1,sig2) + 3.0);
		      dOmega = (value[2] - value[3])/(2.0*coer*hx);
              AngMom = AngMom + 1.0/8.0/(aY[j][i].u*aY[j][i].u)*x*x*dOmega;
              break;
	      }
 	      else if (barrier_rdown1 > 0.9) {
              theta = PetscAtanReal((y + 3.0)/x);
              ierr = indexfunc_down(info, theta, 1.5 + 2.0*coer*hx,
                        1.5 + 0.0*coer*hx, value, aY, user); CHKERRQ(ierr);
		      part = 1.0/(value[1]*PetscExpReal(q_func(x,y,a,sig1,sig2))*1.5);
              //aY[j][i].u = aY[j][i].u;
              aY[j][i].u = value[0] + 4.0/3.0*coer*hx*part*(x*qro(x,y,a,sig1,sig2)
			  	        + (y + 3.0)*qz(x,y,a,sig1,sig2) + 3.0);
              break;
          }
       }
  }
	
  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j + user->D;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hx * i + user->L; 
          if(i == mx - 1 || j == 0 || j == my - 1) {
	      aG[j][i].u = 0.0;
	      aG[j][i].v = 0.0;
	  }
 	  else {
	          uxl = (i==0) ? aY[j][i+1].u : aY[j][i-1].u;
              uxr = (i==mx-1) ? aY[j][i-1].u : aY[j][i+1].u;
              uyl = (j==0) ? aY[j+1][i].u : aY[j-1][i].u;
              uyr = (j==my-1) ? aY[j-1][i].u : aY[j+1][i].u;
              vxl = (i==0) ? aY[j][i+1].v : aY[j][i-1].v;
              vxr = (i==mx-1) ? aY[j][i-1].v : aY[j][i+1].v;
              vyl = (j==0) ? aY[j+1][i].v : aY[j-1][i].v;
              vyr = (j==my-1) ? aY[j-1][i].v : aY[j+1][i].v;

              uxx = (uxl - 2.0 * aY[j][i].u + uxr)/(hx*hx);
              uyy = (uyl - 2.0 * aY[j][i].u + uyr)/(hy*hy);
              u1 = (aY[j][i+1].u - uxl) / (2.0*hx*x);
              u2 = 1.0/4.0*qro2z2(x,y,a,sig1,sig2)*aY[j][i].u;
              u3 = 1.0/16.0*x*x
                     * (PetscPowReal(1.0/(2.0*hx)*(aY[j][i+1].v - vxl),2.0)
                         + PetscPowReal(1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v),2.0))
                     * PetscPowReal(aY[j][i].u, -7.0);

              vxx = (vxl - 2.0 * aY[j][i].v + vxr)/(hx*hx);
              vyy = (vyl - 2.0 * aY[j][i].v + vyr)/(hy*hy);
              v1 = 3.0/x*1.0/(2.0*hx)*(aY[j][i+1].v - vxl);
              v2 = dNdro(x,y)*1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
                        + dNdz(x,y)*1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
              v3 = 6.0/(aY[j][i].u)*1.0/(2.0*hx)*(aY[j][i+1].u - uxl)
                                *1.0/(2.0*hx)*(aY[j][i+1].v - vxl)
                   + 6.0/(aY[j][i].u)*1.0/(2.0*hy)*(aY[j+1][i].u - aY[j-1][i].u)
                                *1.0/(2.0*hy)*(aY[j+1][i].v - aY[j-1][i].v);
	       
              aG[j][i].u = 1.0/1.5*(uxx + uyy + 1.0*(u1 + u3))
                                            *(1.5 - barrier(x, y, 1.5 - coer*hx))
                                + 0.0*barrier(x, y, 1.5 - coer*hx);
              aG[j][i].v = 1.0/1.5*(vxx + vyy + 1.0*v1)
                                            *(1.5 - barrier(x, y, 1.5 - coer*hx))
                                + 0.0*barrier(x, y, 1.5 - coer*hx);
      }
      }
  }
  return 0.; 

}

PetscErrorCode AngMom_ADMMass(DMDALocalInfo *info, DM da, Vec Y, 
                                PetscReal AngMom, AppCtx *user) {
  PetscErrorCode   ierr;
  PetscInt         i, j, n, k = 0;
  PetscReal        x, y, hx, hy, dOmega, barrier_rup1,
                   barrier_thetaup1,theta = 0.0,Mass=0.0,
                   coer = user->coer, value[4],
                   mx = info->mx, my = info->my,
                   a = user->a, sig1 = user->sig1,
                   sig2 = user->sig2,barrier_thetadown1;
  DMDACoor2d       **aC;
  Field            **aY;

  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
  hx = (user->R - user->L) / (PetscReal)(info->mx - 1);
  hy = (user->U - user->D) / (PetscReal)(info->my - 1);
 
  n = (PetscInt)(PETSC_PI/(2.0*coer*hx));
  while (k < n) {
      theta = (k + 1)*2.0*coer*hx - PETSC_PI/2.0;
      for (j = info->ys; j < info->ys + info->ym; j++) {
          y = hy * j + user->D;
          for (i = info->xs; i < info->xs + info->xm; i++) {
	          x = hx * i + user->L;
	          barrier_rup1 = barrier_rdown(x, y, 1.5, coer*hx)
                                *barrier_thetadown(x, y, theta, coer*hx);
              if (barrier_rup1 > 0.6) {
                  ierr = indexfunc_down(info, theta, 1.5 + 2.0*coer*hx,
				        1.5 - 2.0*coer*hx, value, aY, user); CHKERRQ(ierr);
                  dOmega = (value[2] - aY[j][i].v)/(2.0*coer*hx);
                  //AngMom = AngMom + PetscAbsReal(PetscSinReal(theta))
                  //          *2.0*coer*hx; 
                  AngMom = AngMom + 1.0/8.0*1.5*1.5*2.0*coer*hx
                              * PetscPowReal(PetscAbsReal(PetscSinReal(theta)),3.0)
                              * PetscExpReal(q_func(x,y,a,sig1,sig2))*1.5*1.5
                              * dOmega*aY[j][i].u*aY[j][i].u;
                  //            * aY[j][i].u*aY[j][i].u*x*x*dOmega;
                  //PetscPrintf(PETSC_COMM_WORLD,
                //"AngMom=%g, dOmega=%g,theta=%g\n", AngMom, dOmega, theta); CHKERRQ(ierr);
                goto here;
               }
           }    
       }
       here:
       k = k + 1;
  }      
  for (j = info->ys; j < info->ys + info->ym; j++) {
          y = hy * j + user->D;
          for (i = info->xs; i < info->xs + info->xm; i++) {
	          x = hx * i + user->L;
             if (i == mx - 1){
                   Mass = Mass + (aY[j][i].u - aY[j][i-1].u)/hx*hy;
               }
               if (j == 0){
                   Mass = Mass + (aY[j][i].u - aY[j+1][i].u)/hy*hx;
               }
               if (j == my - 1){
                   Mass = Mass + (aY[j][i].u - aY[j-1][i].u)/hy*hx;
               }
 
          }
  }
  PetscPrintf(PETSC_COMM_WORLD,
           "AngMom=%g,Mass=%g\n", AngMom, Mass); CHKERRQ(ierr);

  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}


