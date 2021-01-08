//Codigo Modelo de Evaluacion Anchoveta 2020 version 1
//Autor del codigo por Fernando Espindola
GLOBALS_SECTION
 #include <admodel.h>
 #include <stdio.h>
 #include <time.h>
 #include <limits>
 #include <math.h>
 #include <fvar.hpp>
 time_t start,finish;
 long hour,minute,second;
 double elapsed_time;
// adstring nombre_file;
// adstring nombre_data;
 ofstream mcmc_report("mcmc.csv");

 double eps = std::numeric_limits<double>::epsilon();

 long _stack = 2500000000;

TOP_OF_MAIN_SECTION
 time(&start);
 arrmblsize = 50000000L; 
 gradient_structure::set_NUM_RETURN_ARRAYS(15);
 gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e8); 
 gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e8); 
 gradient_structure::set_MAX_NVAR_OFFSET(50000); 
 gradient_structure::set_NUM_DEPENDENT_VARIABLES(50000); 

 
DATA_SECTION
 int iseed //genera_numero_aleatorio
 !!CLASS random_number_generator rng(iseed);

 init_int ntime
 init_int nedades
 init_int ntallas
 init_matrix matdat(1,ntime,1,18)
 init_vector edades(1,nedades)
 init_vector Tallas(1,ntallas)
 init_matrix Ctot_ch(1,ntime,1,ntallas)
 init_matrix Ctot_pe(1,ntime,1,ntallas)
 init_matrix Ccru(1,ntime,1,ntallas)
 init_matrix Wmed(1,ntime,1,ntallas)
 init_vector msex(1,ntallas)

 !! ad_comm::change_datafile_name("ancho_opt.ctl");

 init_number priorRo
 init_int    fase_Ro

 number log_prior_Ro

 !! log_prior_Ro = log(priorRo);

 init_vector cvar(1,3)
 init_int    fasedevR
 init_int    fasedevNo

 init_vector dt(1,4)

 init_int    nbloques1       // Chilean selectivity
 init_vector ybloques1(1,nbloques1)
 init_int    nbloques2       // Peruvian Selectivity
 init_vector ybloques2(1,nbloques2)
 init_int    nbloques3
 init_vector ybloques3(1,nbloques3)

 init_vector pars_Bio(1,7)
 init_vector cvpars_Bio(1,7) 
 init_int    opt_Linf 
 init_int    opt_k
 init_int    opt_Lo
 init_int    opt_alfa
 init_int    opt_beta
 init_int    opt_M
 init_int    opt_h

 number log_Linfprior
 number log_kprior         
 number log_Loprior        
 number log_alfa_prior
 number log_beta_prior
 number log_M_prior
 number log_h_prior
 number nanos
 number sigmaRo
 number sigmaNo
 number sigmaRt


  !! log_Linfprior= log(pars_Bio(1));
  !! log_kprior = log(pars_Bio(2));
  !! log_Loprior = log(pars_Bio(3));
  !! log_alfa_prior = log(pars_Bio(4));
  !! log_beta_prior = log(pars_Bio(5));
  !! log_M_prior = log(pars_Bio(6));
  !! log_h_prior = log(pars_Bio(7));
  !! nanos = int(ntime/2);
  !! sigmaRo = cvar(1);
  !! sigmaNo = cvar(2);
  !! sigmaRt = cvar(3);

 init_number A50
 init_number cv_A50
 init_number rango
 init_number cv_rango

 init_int    opt_devR
 init_int    opt_Scru1
 init_int    opt1_fase
 init_int    opt_q
 init_int    faseF
 number      logA50_prior
 number      logRango_prior

 !! logA50_prior=log(A50);
 !! logRango_prior=log(rango);

 init_number gama
 init_number opt_Reclu
 init_number cvF
 init_number pbr
 init_int    ntime_sim
 init_int    nrun
 init_int    sem_year
 init_int    sem1
 init_int    sem2
 init_number pR

 number neglog19

 !! neglog19=-1.0*log(19);

 init_int nbloqueqmph
 init_vector ybloqueqmph(1,nbloqueqmph)
 init_int opt_qmph

 init_number low_Len50_S;   init_number upp_Len50_S;   init_int phase_Len50_S;
 init_number low_LenDiff_S; init_number upp_LenDiff_S; init_int phase_LenDiff_S;

 init_vector uno(1,8)
 init_vector dos(1,8)

 number nuno
 number ndos

 !! nuno=double(size_count(uno));
 !! ndos=double(size_count(dos));

// !!cout << msex <<endl;exit(1);
 
 int reporte_mcmc

INITIALIZATION_SECTION
// defino un valor inicial de log_reclutamiento promedio (factor de escala)

 log_Ro              log_prior_Ro
 log_Linf            log_Linfprior
 log_k               log_kprior
 log_Lo              log_Loprior
 log_alfa            log_alfa_prior
 log_beta            log_beta_prior
 log_A50c            logA50_prior
 log_rangoc          logRango_prior
 log_selA            0.09
 log_selB            -1.4
 log_selC            1.5
 log_selD            0.09
 log_selE            -1.4
 log_selF            1.5



PARAMETER_SECTION

 init_bounded_vector log_selA(1,nbloques1,-2.3,1.1,opt1_fase)
 init_bounded_vector log_selB(1,nbloques1,-2.3,1.1,opt1_fase)
 init_bounded_vector log_selC(1,nbloques1,0.1,3.6,opt1_fase)

 init_bounded_vector log_selD(1,nbloques2,-2.3,1.1,opt1_fase)
 init_bounded_vector log_selE(1,nbloques2,-2.3,1.1,opt1_fase)
 init_bounded_vector log_selF(1,nbloques2,0.1,3.6,opt1_fase)


 init_bounded_vector log_A50c(1,nbloques3,-2.3,0.3,opt_Scru1)  
 init_bounded_vector log_rangoc(1,nbloques3,-2.3,-0.2,opt_Scru1)

 init_bounded_vector Len50_S(1,nbloques3,low_Len50_S,upp_Len50_S,phase_Len50_S)
 init_bounded_vector LenDiff_S(1,nbloques3,low_LenDiff_S,upp_LenDiff_S,phase_LenDiff_S)

// Parametros reclutamientos y mortalidades
 init_bounded_number log_Ro(-10,10,fase_Ro)
 init_bounded_dev_vector dev_log_Ro(1,ntime,-10,10,fasedevR)
 init_bounded_dev_vector dev_log_No(1,nedades,-10,10,fasedevNo)
 init_bounded_vector log_Ft_ch(1,ntime,-20,1.39,faseF)
 init_bounded_vector log_Ft_pe(1,ntime,-20,1.39,faseF)

// capturabilidades
 init_number log_qtpe(opt_q)
 init_number log_qtch(opt_q)
 init_vector log_qmph(1,nbloqueqmph,opt_qmph)

// Crecim y M

 init_number log_Linf(opt_Linf)
 init_number log_k(opt_k)
 init_number log_Lo(opt_Lo)
 init_number log_alfa(opt_alfa)
 init_number log_beta(opt_beta)

//Defino las variables de estado 
 vector BMflo(1,ntime);
 vector Bcru1(1,ntime);
 vector Bcru2(1,ntime);
 vector Bmph(1,ntime);
 vector pred_Desemb_ch(1,ntime);
 vector pred_Desemb_pe(1,ntime);
 vector anos(1,ntime);
 vector AcuPE(1,ntime);
 vector AcuCH(1,ntime);
 vector fech(1,ntime);
 vector MPH(1,ntime);
 vector fech_mph(1,ntime);
 vector Desemb_ch(1,ntime);
 vector Desemb_pe(1,ntime);
 vector yCH_fecha(1,ntime);
 vector yPE_fecha(1,ntime);
 vector likeval(1,12);
 vector Neq(1,nedades);
 vector Rpred(1,ntime);
 vector Unos_edad(1,nedades);
 vector Unos_anos(1,ntime);
 vector Unos_tallas(1,ntallas);
 vector mu_edad(1,nedades);
 vector sigma_edad(1,nedades);
 vector BDo(1,ntime);
 vector No(1,nedades)
 vector prior(1,10)
 vector Lmed_Chi_obs(1,ntime)
 vector Lmed_Per_obs(1,ntime)
 vector Lmed_Cru_obs(1,ntime)

 matrix cv_index(1,6,1,ntime)
 matrix nm(1,3,1,ntime)
 matrix Sel_f_ch(1,ntime,1,nedades)
 matrix Sel_f_pe(1,ntime,1,nedades)
 matrix Scru(1,ntime,1,nedades)
 matrix Scru_tallas(1,ntime,1,ntallas)
 matrix Ftot_ch(1,ntime,1,nedades)
 matrix Ftot_pe(1,ntime,1,nedades)
 matrix Ftot(1,ntime,1,nedades)
 matrix Z(1,ntime,1,nedades)
 matrix N(1,ntime,1,nedades)
 matrix NM(1,ntime,1,nedades)
 matrix NMD(1,ntime,1,ntallas)
 matrix NDv(1,ntime,1,ntallas)
 matrix Nrec(1,ntime,1,ntallas)
 matrix NVcru1(1,ntime,1,ntallas)
 matrix NVcru2(1,ntime,1,ntallas)
 matrix NVflo(1,ntime,1,ntallas)
 matrix pred_Ctot_ch(1,ntime,1,ntallas)
 matrix pred_Ctot_pe(1,ntime,1,ntallas)
 matrix Edad_Ctot_ch(1,ntime,1,nedades)
 matrix Edad_Ctot_pe(1,ntime,1,nedades)
 matrix S(1,ntime,1,nedades)
 matrix pobs_fch(1,ntime,1,ntallas)
 matrix ppred_fch(1,ntime,1,ntallas)
 matrix pobs_fpe(1,ntime,1,ntallas)
 matrix ppred_fpe(1,ntime,1,ntallas)
 matrix pobs_cru(1,ntime,1,ntallas)
 matrix ppred_cru(1,ntime,1,ntallas)
 matrix Prob_talla(1,nedades,1,ntallas)
 vector mu_edad1_time(1,ntime)
 vector delta_time(1,ntime)
 4darray MatricesConv(1,5,1,ntime,1,nedades,1,ntallas)
 matrix ClavePro(1,nedades,1,ntallas)
 matrix Nv(1,ntime,1,nedades)
 matrix NMDv(1,ntime,1,nedades)
 matrix NPt(1,ntime,1,ntallas)

 matrix S1(1,nbloques1,1,nedades)
 matrix S2(1,nbloques2,1,nedades)
 matrix S3(1,nbloques3,1,nedades)
 matrix S3_tallas(1,nbloques3,1,ntallas)
 matrix S4(1,nbloques3,1,nedades)

 number suma1
 number suma2
 number suma3
 number suma4
 number suma5
 number suma6
 number res1
 number res2
 number res3
 number res4
 number n1
 number n2
 number n3
 number n4
 number cuenta1
 number cuenta2
 number cuenta3
 number t1
 number t2
 number So
 number Linf
 number linf
 number k
 number ki
 number K
 number Lo
 number M
 number h
 number RPRp
 number phi
 number a
 number sl
 number sr
 number b
 number bl
 number br
 number ide
 number BDp
 number junk

 vector Np(1,nedades)
 vector Zpbr(1,nedades)
 vector Fpbr(1,nedades)
 vector Spbr(1,nedades)
 vector Sp(1,nedades)

 3darray Npy(1,nrun,1,ntime_sim,1,nedades) 
 3darray CTP(1,nrun,1,ntime_sim,1,ntallas)
 3darray NMDpy(1,nrun,1,ntime_sim,1,ntallas)
 matrix BDpy(1,nrun,1,ntime_sim)
 matrix RPRpy(1,nrun,1,ntime_sim)
 matrix YTP(1,nrun,1,ntime_sim)
 vector runo(1,nuno)
 vector rdos(1,ndos)
 number mea1
 number std1
 number min1
 number max1
 number mea2
 number std2
 number min2
 number max2

 objective_function_value f

 sdreport_vector BD(1,ntime) // 
 sdreport_vector BT(1,ntime) // 
 sdreport_vector RPR(1,ntime) // 
 sdreport_vector RPRlp(1,ntime) // 
 sdreport_vector Reclutas(1,ntime) // 
// sdreport_number CBA 

 sdreport_number SSBo
 sdreport_number alfa
 sdreport_number beta

 sdreport_vector Ftot_ref(1,ntime)
 sdreport_vector BD_anual(1,nanos)
 sdreport_vector FT_anual(1,nanos)
 sdreport_vector Lmed_Chi_pre(1,ntime)
 sdreport_vector Lmed_Per_pre(1,ntime)
 sdreport_vector Lmed_Cru_pre(1,ntime)
 sdreport_vector pred_AcuPE(1,ntime)
 sdreport_vector pred_AcuCH(1,ntime)
 sdreport_vector pred_MPH(1,ntime)
 
 
 likeprof_number F_fin;
 likeprof_number Deple_fin;
 
 vector log_Fch_prior(1,ntime)
 vector log_Fpe_prior(1,ntime)
 

PRELIMINARY_CALCS_SECTION
// leo la matriz de indices

 anos=column(matdat,1);
 AcuPE=column(matdat,2);
 cv_index(1)=column(matdat,3);
 AcuCH=column(matdat,4);
 cv_index(2)=column(matdat,5);
 fech=column(matdat,6);
 MPH=column(matdat,7);
 cv_index(4)=column(matdat,8);
 Desemb_ch=column(matdat,9);
 cv_index(5)=column(matdat,10);
 Desemb_pe=column(matdat,11);
 cv_index(6)=column(matdat,12);
 nm(1)=column(matdat,13);
 nm(2)=column(matdat,14);
 nm(3)=column(matdat,15);
 yCH_fecha=column(matdat,16);
 yPE_fecha=column(matdat,17);
 fech_mph=column(matdat,18);


// junk=trun_rnorm(5,4,2,8);
// cout << "MEAN: " << junk << endl;exit(1);

//   cout << "log_Ro=" << log_Ro << endl;exit(1);
//   cout << "VF_year:" << VF_year << endl;exit(1);

 Unos_edad=1;// lo uso en  operaciones matriciales con age
 Unos_anos=1;// lo uso en operaciones matriciales con year
 Unos_tallas=1;// lo uso en operaciones matriciales con year
 reporte_mcmc=1;

RUNTIME_SECTION

  convergence_criteria 1.e-01,1.e-02,1.e-03,1e-5,1e-6
  maximum_function_evaluations 100,100,200,3000,5000


PROCEDURE_SECTION
// se listan las funciones que contienen los calculos

 Eval_prob_talla_edad_variable();
 Eval_selectividad_normal();
 Eval_mortalidades();
 Eval_abundancia();
 Eval_biomasas();
 Eval_capturas_predichas();
 Eval_indices();
 Eval_deinteres();
 Eval_logverosim();
 Eval_funcion_objetivo();
 Eval_mcmc();
// Eval_CTP();

 if (last_phase()) {Eval_CTP();}


FUNCTION Eval_prob_talla_edad_variable
 int i,j,m;

 //aqui calcula la Matriz de Conversion variable en el tiempo
 Linf=mfexp(log_Linfprior);
 k=mfexp(log_kprior);
 Lo=mfexp(log_Lo);
 //cout << "k = " << k << endl;exit(1);

  for(m=1;m<=5;m++){ 
   switch (m) 
    { 
       case 1: delta_time= yCH_fecha; //Chile capturas
                break; 
       case 2: delta_time= yPE_fecha; //Peru Capturas
                break; 
       case 3: delta_time= fech; // Acustica Chile
               break; 
       case 4: delta_time=dt(4); // Acustica Peru
               break;
       case 5: delta_time=fech_mph;// MPH
               break;
       default: printf("ALK no variable");
                delta_time=0.0; 
                break;   
    }
    if(m==1 || m==2 || m==3 || m==5){
     mu_edad1_time=Lo+(Linf-Lo)*(1.0-mfexp(-k*delta_time));
      for(j=1;j<=ntime;j++){
        mu_edad(1)=mu_edad1_time(j);
        for(i=2;i<=nedades;i++){
          mu_edad(i)=mu_edad(i-1)*mfexp(-k)+Linf*(1-mfexp(-k));
        }
        sigma_edad=mfexp(log_alfa)+mfexp(log_beta)*mu_edad;
        MatricesConv(m,j)=ALK(mu_edad,sigma_edad,Tallas);
      }
    } else {
     mu_edad(1)=Lo+(Linf-Lo)*(1.0-mfexp(-k*dt(4)));
     for(j=1;j<=ntime;j++){
      for(i=2;i<=nedades;i++){
       mu_edad(i)=mu_edad(i-1)*mfexp(-k)+Linf*(1-mfexp(-k));
       sigma_edad=mfexp(log_alfa)+mfexp(log_beta)*mu_edad;
       MatricesConv(4,j)=ALK(mu_edad,sigma_edad,Tallas);
      }
     }
    }
  }
//VB tradicional
  mu_edad(1)=Lo;
  for(int i=2;i<=nedades;i++){
    mu_edad(i)=mu_edad(i-1)*mfexp(-k)+Linf*(1-mfexp(-k));
  }
  sigma_edad=mfexp(log_alfa)+mfexp(log_beta)*mu_edad;
  Prob_talla=ALK(mu_edad,sigma_edad,Tallas);


// cout << "Prob_talla" << Prob_talla(1) << endl;exit(1);
//-------------------------------------------------------
FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
        int i, j;
        dvariable z1;
        dvariable z2;
        int si,ni; si=mu.indexmin(); ni=mu.indexmax();
        int sj,nj; sj=x.indexmin(); nj=x.indexmax();
        dvar_matrix pdf(si,ni,sj,nj);
        pdf.initialize();
        double xs=0.5*(x[sj+1]-x[sj]);
        for(i=si;i<=ni;i++){
          for(j=sj;j<=nj;j++){
            z1=((x(j)-xs)-mu(i))/sig(i);
            z2=((x(j)+xs)-mu(i))/sig(i);
            pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);

          }
          pdf(i)/=sum(pdf(i));
        }
        return(pdf);
//-------------------------------------------------------

// cout << "Prob_talla" << Prob_talla << endl;exit(1);
// cout << "mu_edad" << mu_edad << endl;exit(1);
// cout << "sigma_edad" << sigma_edad << endl;exit(1); 
// cout << "Tallas" << Tallas << endl;exit(1);


FUNCTION Eval_selectividad_normal
 int i,j;

 for(j=1;j<=nbloques1;j++){
    a  = mfexp(log_selA(j));
    sl = mfexp(log_selB(j));
    sr = mfexp(log_selC(j));
    for(i=1;i<=nedades;i++){
       if(i<=a){S1(j,i)=pow(2,-1*square((i-a)/sl));}
       else{S1(j,i)=pow(2,-1*square((i-a)/sr));}
    }
 S1(j)=S1(j)/(max(S1(j))+eps);
 }

 for(j=1;j<=nbloques2;j++){
    b  = mfexp(log_selD(j));
    bl = mfexp(log_selE(j));
    br = mfexp(log_selF(j));
    for(i=1;i<=nedades;i++){
       if(i<=b){S2(j,i)=pow(2,-1*square((i-b)/bl));}
       else{S2(j,i)=pow(2,-1*square((i-b)/br));}
    }
 S2(j)=S2(j)/(max(S2(j))+eps);
 }

 for(i=1;i<=ntime;i++){
   for(j=1;j<=nbloques1;j++){if(anos(i)>=ybloques1(j)){Sel_f_ch(i)=S1(j);}}
 }
 for(i=1;i<=ntime;i++){
   for(j=1;j<=nbloques2;j++){if(anos(i)>=ybloques2(j)){Sel_f_pe(i)=S2(j);}}
 }

// Selectividad de cruceros
 for (j=1;j<=nbloques3;j++){
  S3(j)=1./(1.+mfexp(neglog19*(edades-mfexp(log_A50c(j)))/mfexp(log_rangoc(j))));
 }
 for(i=1;i<=ntime;i++){
   for(j=1;j<=nbloques3;j++){if (anos(i)>=ybloques3(j)){Scru(i)=S3(j);}}
 }

//Selectividad de cruceros acusticos a la talla
 for(j=1;j<=nbloques3;j++)
  S3_tallas(j)=1./(1.+mfexp(neglog19*(Tallas-mfexp(Len50_S(j)))/mfexp(LenDiff_S(j))));
 for(int i=1;i<=ntime;i++){
  for(int j=1;j<=nbloques3;j++){if (anos(i)>=ybloques3(j)){Scru_tallas(i)=S3_tallas(j);}}
 }

// cout << "log_rangoc" << log_rangoc << endl;exit(1);

FUNCTION Eval_mortalidades

 M=mfexp(log_M_prior);

 Ftot_ch=elem_prod(Sel_f_ch,outer_prod(mfexp(log_Ft_ch),Unos_edad));
 Ftot_pe=elem_prod(Sel_f_pe,outer_prod(mfexp(log_Ft_pe),Unos_edad));
 Ftot=Ftot_ch+Ftot_pe;

 Z=Ftot+M;
 S=mfexp(-1.0*Z);

 // F promedio dos ultimso penultimo semestres
 
FUNCTION Eval_abundancia
 int i, j;

 // Biomasa desovante virgen de largo plazo
 No(1)=mfexp(log_Ro+0.5*square(sigmaRo));
 for (int j=2;j<=nedades;j++){
   No(j)=No(j-1)*mfexp(-1.*M);
 }
 No(nedades)=No(nedades)/(1-mfexp(-1.*M));

// No(nedades)=No(nedades)*mfexp(-1.*M)/(1-mfexp(-1.*M));
// cout << "No" << No << endl;exit(1);

 SSBo=sum(elem_prod(No*mfexp(-dt(4)*M)*Prob_talla,elem_prod(msex,Wmed(1))));
 phi=SSBo/mfexp(log_Ro);
 // Stock-recluta
 h=mfexp(log_h_prior);
//
 if(opt_Reclu==1 || opt_Reclu==4){//Beverton-Holt
  alfa = 4*h*mfexp(log_Ro)/(5*h-1);
  beta = (1-h)*SSBo/(5*h-1);//  
 } else if(opt_Reclu==2){//Ricker
  alfa = 1.25*log(5*h)-log(phi) ; 
  beta = 1.25*log(5*h)/SSBo;
 }

//Abundancia inicial de equilibrio como referencia para el primer year
// Neq(1)=mfexp(log_Ro)*No(1);

 Neq(1)=mfexp(log_Ro+0.5*square(sigmaRo));
 for (int j=2;j<=nedades;j++){
   Neq(j)=Neq(j-1)*mfexp(-Z(1,j-1));   
 }
 Neq(nedades)=Neq(nedades)/(1-mfexp(-1.*Z(1,nedades)));

// Neq(nedades)=Neq(nedades)*mfexp(-1.*Z(1,nedades))/(1-mfexp(-1.*Z(1,nedades)));
 N(1)=elem_prod(Neq,mfexp(dev_log_No));

 if(opt_Reclu==1 || opt_Reclu==2){
  BD(1)=sum(elem_prod(elem_prod(N(1),mfexp(-dt(4)*Z(1)))*Prob_talla,elem_prod(Wmed(1),msex)));
  Rpred(1)=N(1,1);
 } else if(opt_Reclu==3){
 for(int i=2;i<=ntime;i++){//de la edad 2 en adelante 
   N(i,1)=mfexp(log_Ro+dev_log_Ro(i)+0.5*square(sigmaRo));
  }
  BD(1)=sum(elem_prod(elem_prod(N(1),mfexp(-dt(4)*Z(1)))*Prob_talla,elem_prod(Wmed(1),msex)));
  Rpred(1)=N(1,1);
 }

 for(int i=1;i<ntime;i++){

  if(opt_Reclu==1){//Opcion_Beverton_Holt_Model
    Rpred(i+1)=alfa*BD(i)/(beta+BD(i));
    N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i+1));
  }
  else if(opt_Reclu==2){//Opcion_Ricker_Model
    Rpred(i+1)=BD(i)*mfexp(alfa-1.0*beta*BD(i));
//    Rpred(i+1)=(log_Ro*BD(i)/SSBo)*mfexp(log(5*h)*((1-BD(i)/SSBo)/0.8));
    N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i+1));
  }
  else if(opt_Reclu==3){
    N(i,1)=N(i,1);
  }
  else if(opt_Reclu==4){
    Rpred(i+1)=alfa*BD(i)/(beta+pow(BD(i),gama));
    N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i+1));
  }
  //N(i+1,1)=Rpred(i+1)*mfexp(dev_log_Ro(i+1));
  N(i+1)(2,nedades)=++elem_prod(N(i)(1,nedades-1),S(i)(1,nedades-1));
  N(i+1,nedades)=N(i+1,nedades)+N(i,nedades)*S(i,nedades);// grupo plus
  BD(i+1)=sum(elem_prod(elem_prod(N(i+1),mfexp(-dt(4)*Z(i+1)))*Prob_talla,elem_prod(Wmed(i+1),msex)));
 }
// cout << "msex:"  << msex << endl;exit(1);

FUNCTION Eval_biomasas

// matrices y vectores de abundancias derivadas
 Reclutas=column(N,1);

// NVcru1=elem_prod(elem_prod(N,mfexp(elem_prod(outer_prod(fech,Unos_edad),Z))),Scru)*Prob_talla;//acoustic chile
 for(int j=1;j<=ntime;j++){
  NVcru1(j)=elem_prod(elem_prod(N(j),mfexp(fech(j)*Z(j))),Scru(j))*MatricesConv(3,j);//acoustic chile
 }
 NVcru1=elem_prod(NVcru1,Scru_tallas);

 NVcru2=elem_prod(elem_prod(N,mfexp(-dt(1)*Z)),Scru)*MatricesConv(4,1);//acoustic peru

// cout << "fech:" << fech << endl;exit(1);

// NMD=elem_prod(N,mfexp(-dt(4)*Z))*Prob_talla;//spawning season
// NMD=elem_prod(NMD,outer_prod(Unos_anos,msex));
 for(int j=1;j<=ntime;j++){
  NMD(j)=elem_prod(N(j),mfexp(fech_mph(j)*Z(j)))*MatricesConv(5,j);//spawning season
 }
 NMD=elem_prod(NMD,outer_prod(Unos_anos,msex));

 for(int j=1;j<=ntime;j++){
  Bmph(j)=NMD(j)*Wmed(j);
 }

// vectores de biomasas derivadas
 NPt=N*Prob_talla;

 for(int j=1;j<=ntime;j++){
  Bcru1(j)=NVcru1(j)*Wmed(j);
  Bcru2(j)=NVcru2(j)*Wmed(j);
  BT(j)=NPt(j)*Wmed(j);  
 }


FUNCTION Eval_capturas_predichas

// matrices de capturas predichas por edad y year
 Edad_Ctot_ch=elem_prod(elem_div(Ftot_ch,Z),elem_prod(1.-S,N));
 Edad_Ctot_pe=elem_prod(elem_div(Ftot_pe,Z),elem_prod(1.-S,N));
 for (int i=1;i<=ntime;i++){pred_Ctot_ch(i)=Edad_Ctot_ch(i)*MatricesConv(1,i);}
 for (int i=1;i<=ntime;i++){pred_Ctot_pe(i)=Edad_Ctot_pe(i)*MatricesConv(2,i);}

// vectores de desembarques predichos por year
 for(int j=1;j<=ntime;j++){
  pred_Desemb_ch(j)=pred_Ctot_ch(j)*Wmed(j);
  pred_Desemb_pe(j)=pred_Ctot_pe(j)*Wmed(j);
 }

// matrices de proporcion de capturas por edad y year
 pobs_fch=elem_div(Ctot_ch,outer_prod(rowsum(Ctot_ch+eps),Unos_tallas));
 ppred_fch=elem_div(pred_Ctot_ch,outer_prod(rowsum(pred_Ctot_ch+eps),Unos_tallas));

 pobs_fpe=elem_div(Ctot_pe,outer_prod(rowsum(Ctot_pe+eps),Unos_tallas));
 ppred_fpe=elem_div(pred_Ctot_pe,outer_prod(rowsum(pred_Ctot_pe+eps),Unos_tallas));

 pobs_cru=elem_div(Ccru,outer_prod(rowsum(Ccru+eps),Unos_tallas));
 ppred_cru=elem_div(NVcru1,outer_prod(rowsum(NVcru1),Unos_tallas));


FUNCTION Eval_indices
 int i,j;

 pred_AcuPE=mfexp(log_qtpe)*Bcru2;
 pred_AcuCH=mfexp(log_qtch)*Bcru1;
// pred_MPH=mfexp(log_qmph)*BD;
 for (int i=1;i<=ntime;i++){
   for (int j=1;j<=nbloqueqmph;j++){
      if (anos(i)>=ybloqueqmph(j)){
        pred_MPH(i)=mfexp(log_qmph(j))*Bmph(i);}
   }
 }



FUNCTION Eval_deinteres
// aca pongo las variables de interes para mcmc y otras yerbas

// Rutina para calcular RPR
 Nv=N;// solo para empezar los calculos
  
 for (int i=1;i<ntime;i++){
   Nv(i+1)(2,nedades)=++Nv(i)(1,nedades-1)*mfexp(-1.0*M);
 }
 for (int i=2;i<ntime;i++){
   Nv(i,nedades)=Nv(i,nedades)+mfexp(-1.0*M)*Nv(i,nedades);// grupo plus
 }

 NDv=elem_prod((Nv*mfexp(-dt(4)*M))*Prob_talla,outer_prod(Unos_anos,msex));
 for(int j=1;j<=ntime;j++){
  BDo(j)=NDv(j)*Wmed(j);
 }
 RPR=elem_div(BD,BDo);
 RPRlp=BD/SSBo;
 Deple_fin=0.5*(RPRlp(ntime-1)+RPRlp(ntime));// F promedio dos ultimso penultimo semestres
 

FUNCTION Eval_logverosim
// esta funcion evalua el nucleo de las -log-verosimilitudes marginales para
// series con datos 0.
 int i;

 suma1=0;  suma2=0;  suma4=0;  suma5=0;  suma6=0;
 res1=0; res2=0; res3=0; res4=0;
 n1=0; n2=0; n3=0; n4=0;

 for (i=1;i<=ntime;i++){
  if (AcuPE(i)>0){
   suma1+=square((log(AcuPE(i))-log(pred_AcuPE(i)))/cv_index(1,i));
   res1+=square(log(AcuPE(i))-log(pred_AcuPE(i)));
   n1+=1;}
  if (AcuCH(i)>0){
   suma2+=square((log(AcuCH(i))-log(pred_AcuCH(i)))/cv_index(2,i));
   res2+=square(log(AcuCH(i))-log(pred_AcuCH(i)));
   n2+=1;}
  if (MPH(i)>0){
   suma4+=square((log(MPH(i))-log(pred_MPH(i)))/cv_index(4,i));
   res4+=square(log(MPH(i))-log(pred_MPH(i)));
   n4+=1;}
  if (Desemb_ch(i)>0){
   suma5+=square((log(Desemb_ch(i)+eps)-log(pred_Desemb_ch(i)))/cv_index(5,i));}
  if (Desemb_pe(i)>0){
   suma6+=square((log(Desemb_pe(i)+eps)-log(pred_Desemb_pe(i)))/cv_index(6,i));}
 }


FUNCTION Eval_funcion_objetivo
// se calcula la F.O. como la suma de las -logver
// lognormal

 likeval(1)=0.5*suma1;//acustica_Peru
 likeval(2)=0.5*suma2;//acustica_Chile
// likeval(3)=0.5*suma3;//rpe
 likeval(4)=0.5*suma4;//MPH
 likeval(5)=0.5*suma5;//Desemb_ch
 likeval(6)=0.5*suma6;//Desemb_pe

// multinomial

 likeval(7)=-1.0*sum(elem_prod(elem_prod(pobs_fch,outer_prod(nm(1),Unos_tallas)),log(ppred_fch+eps)));
 likeval(8)=-1.0*sum(elem_prod(elem_prod(pobs_fpe,outer_prod(nm(2),Unos_tallas)),log(ppred_fpe+eps)));
 likeval(9)=-1.0*sum(elem_prod(elem_prod(pobs_cru,outer_prod(nm(3),Unos_tallas)),log(ppred_cru+eps)));

// lognormal Ninicial y Reclutas y Proporcion reclutas
 likeval(10)=1./(2*square(sigmaRo))*norm2(dev_log_Ro);
 likeval(11)=1./(2*square(sigmaNo))*norm2(dev_log_No);
 
// prioris
 prior(1)=1./(2*square(cvpars_Bio(1)))*square(log_Linf-log_Linfprior);
 prior(2)=1./(2*square(cvpars_Bio(2)))*square(log_k-log_kprior);
 prior(3)=1./(2*square(cvpars_Bio(3)))*square(log_Lo-log_Loprior);
 prior(4)=1./(2*square(cv_A50))*norm2(log_A50c-logA50_prior);
 prior(5)=1./(2*square(cv_rango))*norm2(log_rangoc-logRango_prior);
// prior(6)=1./(2*square(cv_rango))*norm2(log_rangof_ch-logRango_prior);
// prior(7)=1./(2*square(cv_rango))*norm2(log_rangof_pe-logRango_prior);

  f=sum(likeval)+sum(prior);

 
 if(last_phase){

  Lmed_Chi_obs=Tallas*trans(pobs_fch);
  Lmed_Per_obs=Tallas*trans(pobs_fpe);
  Lmed_Cru_obs=Tallas*trans(pobs_cru);
  Lmed_Chi_pre=Tallas*trans(ppred_fch);
  Lmed_Per_pre=Tallas*trans(ppred_fpe);
  Lmed_Cru_pre=Tallas*trans(ppred_cru);

  for(int i=1;i<=ntime;i++){
   Ftot_ref(i)=mfexp(log_Ft_ch(i))+mfexp(log_Ft_pe(i));
  }
  double cont;
  cont=1;
  for(int j=1;j<ntime;j+=2){
   BD_anual(cont)=(BD(j)+BD(j+1))/2;
   FT_anual(cont)=(Ftot_ref(j)+Ftot_ref(j+1))/2;
   cont=cont+1;
  }

 }
 F_fin=0.5*(Ftot_ref(ntime-1)+Ftot_ref(ntime));



FUNCTION Eval_CTP
 
 for(int r=1;r<=nuno;r++){runo(r)=Reclutas(uno(r));}
 mea1=mean(runo);std1=std_dev(runo);min1=min(runo);max1=max(runo);
 for(int s=1;s<=ndos;s++){rdos(s)=Reclutas(dos(s));}
 mea2=mean(rdos);std2=std_dev(rdos);min2=min(rdos);max2=max(rdos);
// cout << "std1: " << std1 << endl;exit(1);
 //for(int i=1;i<=ntime_sim;i++){
//  for(int j=1;j<=ntallas,j++){

 // }
 //}


 int j=1;
 while(j<=nrun){
   Np=N(ntime);
   Sp=S(ntime);
   for(int i=1;i<=ntime_sim;i++){//semestre_++
     if(i==1){
      Npy(j)(i)(2,nedades)=++elem_prod(Np(1,nedades-1),Sp(1,nedades-1));
      Npy(j)(i)(nedades)+=Np(nedades)*Sp(nedades);
     } else {
      Npy(j)(i)(2,nedades)=++elem_prod(Npy(j)(i-1)(1,nedades-1),Sp(1,nedades-1));
      Npy(j)(i)(nedades)+=Npy(j)(i-1)(nedades)*Sp(nedades);
     }
// cout << "NPy(1,1): " << Npy(1,1) << endl;exit(1);
     if(i%2==0){
      do_trun_rnorm(mea2,std2,min2,max2,j,i);
     }
     if(i%2!=0){
      do_trun_rnorm(mea1,std1,min1,max1,j,i);
     }
     Fpbr=(Ftot(ntime)/max(Ftot(ntime)))*pbr;
     Zpbr=Fpbr+M;
     Sp=exp(-1.*Zpbr);
     NMDpy(j)(i)(1,ntallas)=elem_prod(Npy(j)(i),mfexp(-dt(4)*Zpbr))*MatricesConv(1,ntime);//Prob_talla;
     BDpy(j)(i)=sum(elem_prod(elem_prod(NMDpy(j)(i),msex),Wmed(ntime)));
     CTP(j)(i)(1,ntallas)=elem_prod(elem_div(Fpbr,Zpbr),elem_prod(Npy(j)(i),(1-Spbr)))*MatricesConv(1,ntime);//*Prob_talla;
     YTP(j)(i)=sum(elem_prod(CTP(j)(i),Wmed(ntime)));
     RPRpy(j)(i)=BDpy(j)(i)/SSBo;
   }
 j=j+1;
 }
//
// cout << "Npy: " << Npy(1) << endl;exit(1);
//
//---------------------------------------------------------------
FUNCTION void do_trun_rnorm(dvariable ma, dvariable std, dvariable min, dvariable max, int j, int i)

  dvariable trn;
  trn=ma+randn(rng)*std;
  while(trn<min || trn>max){
   trn=ma+randn(rng)*std;
  }
  Npy(j)(i)(1)=pR*trn;
//----------------------------------------------------------------


REPORT_SECTION

 report << "Year" << endl;
 report << anos << endl;
 report << "AcuPE_obs" << endl;
 report << AcuPE << endl;
 report << "AcuPE_pre" << endl;
 report << pred_AcuPE << endl;
 report << "cv_acus" << endl;
 report << cv_index(1) << endl;
 report << "AcuCH_obs" << endl;
 report << AcuCH << endl;
 report << "AcuCH_pre" << endl;
 report << pred_AcuCH << endl;
 report << "cv_AcuCH" << endl;
 report << cv_index(2) << endl;
 report << "MPH_obs" << endl;
 report << MPH << endl;
 report << "MPH_pre" << endl;
 report << pred_MPH << endl;
 report << "cv_MPH" << endl;
 report << cv_index(4) << endl;
 report << "Desem_Chi_obs" << endl;
 report << Desemb_ch << endl;
 report << "Desem_Chi_pre" << endl;
 report << pred_Desemb_ch << endl;
 report << "cv_Desemb_ch" << endl;
 report << cv_index(5) << endl;
 report << "Desem_Per_obs" << endl;
 report << Desemb_pe << endl;
 report << "Desem_Per_pre" << endl;
 report << pred_Desemb_pe << endl;
 report << "cv_Desemb_pe" << endl;
 report << cv_index(6) << endl;
 report << "Ctot_Chi_obs" << endl;
 report << Ctot_ch << endl;
 report << "Ctot_Per_obs" << endl;
 report << Ctot_pe << endl;
 report << "Ctot_Chi_pre" << endl;
 report << pred_Ctot_ch << endl;
 report << "Ctot_Per_pre" << endl;
 report << pred_Ctot_pe << endl;
 report << "Lmed_Chi_obs" << endl;
 report << Lmed_Chi_obs << endl;
 report << "Lmed_Chi_pre" << endl;
 report << Lmed_Chi_pre << endl;
 report << "Lmed_Per_obs" << endl;
 report << Lmed_Per_obs << endl;
 report << "Lmed_Per_pre" << endl;
 report << Lmed_Per_pre << endl;
 report << "Lmed_Cru_obs" << endl;
 report << Lmed_Cru_obs << endl;
 report << "Lmed_Cru_pre" << endl;
 report << Lmed_Cru_pre << endl;
 report << "Biom_deso" << endl;
 report << BD << endl;
 report << "Biom_cero" << endl;
 report << BDo << endl;
 report << "Biom_total" << endl;
 report << BT << endl;
 report << "SSBo" << endl;
 report << SSBo << endl;
 report << "log_Ro" << endl;
 report << exp(log_Ro) << endl;
 report << "PHI" << endl;
 report << phi << endl;
 report << "Reclu_pre" << endl;
 report << Rpred << endl;
 report << "Reclu_obs" << endl;
 report << Reclutas << endl;
 report << "Reclu_dev" << endl;
 report << dev_log_Ro << endl;
 report << "No_dev" << endl;
 report << dev_log_No << endl;
 report << "Opcion_Reclu" << endl;
 report << opt_Reclu << endl;
 report << "RPR_dina" << endl;
 report << RPR << endl;
 report << "RPR_long" << endl;
 report << RPRlp << endl;
 report << "F_Chile" << endl;
 report << exp(log_Ft_ch) << endl;
 report << "F_Peru" << endl;
 report << exp(log_Ft_pe) << endl;
 report << "F_Total" << endl;
 report <<  Ftot_ref << endl;
 report << "BD_anual" << endl;
 report << BD_anual << endl;
 report << "FT_anual" << endl;
 report << FT_anual << endl;
 report << "F_fin" << endl;
 report << F_fin << endl;
 report << "Ftot" << endl;
 report << Ftot << endl;
 report << "Selec_flo_Chi" << endl;
 report << Sel_f_ch << endl;
 report << "Selec_flo_Per" << endl;
 report << Sel_f_pe << endl;
 report << "Selec_cru" << endl;
 report << Scru << endl;
 report << "Edades" << endl;
 report << edades << endl;
 report << "F_ch" << endl;
 report << Ftot_ch << endl;
 report << "F_pe" << endl;
 report << Ftot_pe << endl;
 report << "Pro_flo_chi_obs" << endl;
 report << pobs_fch << endl;
 report << "Pro_flo_chi_pre" << endl;
 report << ppred_fch << endl;
 report << "sz_pro_flo_chi" << endl;
 report << nm(1) << endl;
 report << "Pro_flo_per_obs" << endl;
 report << pobs_fpe << endl;
 report << "Pro_flo_per_pre" << endl;
 report << ppred_fpe << endl;
 report << "sz_pro_flo_per" << endl;
 report << nm(2) << endl;
 report << "Pro_cru_obs" << endl;
 report << pobs_cru << endl;
 report << "Pro_cru_pre" << endl;
 report << ppred_cru << endl;
 report << "Edad_Ctot_ch" << endl;
 report << Edad_Ctot_ch << endl;
 report << "Edad_Ctot_pe" << endl;
 report << Edad_Ctot_pe << endl;
 report << "sz_pro_cru_obs" << endl;
 report << nm(3) << endl;
 report << "Tallas" << endl;
 report << Tallas << endl;
 report << "Sobre" << endl;
 report << S << endl;
 report << "N" << endl;
 report << N << endl;
 report << "NVcru1" << endl;
 report << NVcru1 << endl;
 report << "Neq" << endl;
 report << Neq << endl;
 report << "No" << endl;
 report << No << endl;
 report << "log_like" << endl;
 report << likeval << endl;
 report << "log_prior" << endl;
 report << prior << endl;
 report << "Prob_talla" << endl;
 report << Prob_talla << endl;
 report << "MatricesConv" << endl;
 report << MatricesConv << endl;
 report << "ClaveCH" << endl;
 report << MatricesConv(1,ntime) << endl;
 report << "ClavePE" << endl;
 report << MatricesConv(2,ntime) << endl;
 report << "ClaveCr" << endl;
 report << MatricesConv(3,ntime) << endl;
 report << "Mu_edad" << endl;
 report << mu_edad << endl;
 report << "Sigma_edad" << endl;
 report << sigma_edad << endl;
 report << "Linf" << endl;
 report << Linf << endl;
 report << "k" << endl;
 report << k << endl;
 report << "Lo" << endl;
 report << Lo << endl;
 report << "alfa" << endl;
 report << alfa << endl;
 report << "beta" << endl;
 report << beta << endl;
 report << "M" << endl;
 report << M << endl;
 report << "h" << endl;
 report << h << endl;
 report << "Peso_medio" << endl;
 report << Wmed << endl;
 report << "Madu_talla" << endl;
 report << msex << endl;
 report << "Madu_edad" << endl;
 report << msex*trans(Prob_talla) << endl;
 report << "Peso_medio_edad" << endl;
 report << Wmed*trans(Prob_talla) << endl;
 report << "AcuPE_cv" << endl;
 report << sqrt(1/n1*res1) << endl;
 report << "AcuCH_cv" << endl;
 report << sqrt(1/n2*res2) << endl;
 report << "MPH_cv" << endl;
 report << sqrt(1/n4*res4) << endl;
//ESTIMA nm 
//
  suma1=0; suma2=0;n1=1;n2=1;n3=1;cuenta1=0;cuenta2=0;cuenta3=0;

  for (int i=1;i<=ntime;i++){ //
   if (sum(pobs_fch(i))>0){
      suma1=sum(elem_prod(ppred_fch(i),1-ppred_fch(i)));
      suma2=norm2(pobs_fch(i)-ppred_fch(i));
      n1=n1*suma1/suma2;
      cuenta1+=1;
   }
   if (sum(pobs_fpe(i))>0){
      suma1=sum(elem_prod(ppred_fpe(i),1-ppred_fpe(i)));
      suma2=norm2(pobs_fpe(i)-ppred_fpe(i));
      n2=n2*suma1/suma2;
      cuenta2+=1;
   }
   if (sum(pobs_cru(i))>0){
      suma1=sum(elem_prod(ppred_cru(i),1-ppred_cru(i)));
      suma2=norm2(pobs_cru(i)-ppred_cru(i));
      n3=n3*suma1/suma2;
      cuenta3+=1;
   }
  }
//
//FIN nm
 report << "Tama_post_flo_Chi" << endl;
 report << pow(n1,1/cuenta1) << endl;
 report << "Tama_post_flo_Per" << endl;
 report << pow(n2,1/cuenta2) << endl;
 report << "Tama_post_Cru" << endl;
 report << pow(n3,1/cuenta3) << endl;
 report << "alfa_edad" << endl;
 report << exp(log_alfa) << endl;
 report << "beta_edad" << endl;
 report << exp(log_beta) << endl;
 report << "post_A50cru" << endl;
 report << exp(log_A50c) << endl;
//ESTIMA CTP
//
//FIN CTP
 report << "PBR" << endl;
 report << pbr << endl;
 report << "nyear_sim" << endl;
 report << ntime_sim << endl;
 report << "nrun_sim" << endl;
 report << nrun << endl;
 report << "Np" << endl;
 report << Np << endl;
 report << "Nproy" << endl;
 report << Npy << endl;
 report << "CTP_proy" << endl;
 report << CTP << endl;
 report << "YTP_proy" << endl;
 report << YTP << endl;
 report << "SSB_proy" << endl;
 report << BDpy << endl; 
 report << "RPR_proy" << endl;
 report << RPRpy << endl;
 

FUNCTION Eval_mcmc

 if(reporte_mcmc == 0)

 mcmc_report << "BD Ft" << endl;//imprime el encabezado
 mcmc_report << BD(ntime) << "," << Ftot_ref(ntime) << endl;//imprime vars

 reporte_mcmc++;


FINAL_SECTION

 time(&finish);
 elapsed_time=difftime(finish,start);
 hour=long(elapsed_time)/3600;
 minute=long(elapsed_time)%3600/60;
 second=(long(elapsed_time)%3600)%60;
 cout<<endl<<endl<<"*********************************************"<<endl;
 cout<<"--Start time:  "<<ctime(&start)<<endl;
 cout<<"--Finish time: "<<ctime(&finish)<<endl;
 cout<<"--Runtime: ";
 cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
 cout<<" "<<endl;
 cout<<"*********************************************"<<endl;
 cout<<" Modelo de Evaluacion de Stock Anchoveta"<<endl;
 cout<<"   del Sur de Peru y Norte de Chile    "<<endl;
 cout<<"                 por                    "<<endl;
 cout<<"          Fernando Espindola            "<<endl;
 cout<<"     INSTITUTO DE FOMENTO PESQUERO      "<<endl;
 cout<<"*********************************************"<<endl;


