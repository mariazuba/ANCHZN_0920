#__________________________________________________________________
#Ro
187.8793
#fase_Ro
1
#cv_recruitments_and_intial_condition_(No)
#sigmaRo(0.62)_y_sigmaNo(0.72)_y_sigmaRt_[sigmaR=0.62_for_B-H_&_Ricker,0.45_for_Desvios]
0.64	0.65	0.3
#fases_para_dev_R_y_dev_No
2
2
#Dt_de_cada_crucero_(mes/6)_es_semestral
#AcuPE(0.83)__AcuCH___RecPE___SSB_(septiembre)
0.83	0.0	0.0	0.16
#______________________________________________________________________
#No_bloques_select_Flo_Chil_anos_de_inicio(primero_siempre_es_el_ano_de_partida)
6
1986	1991	1998.5	2000	2010	2014.5
#No_bloques_select_Flo_Peru_anos_de_inicio(primero_siempre_es_el_ano_de_partida)
7
1986	1989.5	1993.5	1996.5	1998.5	2004	2013
#No_bloques_select_Acustico_(primero_siempre_es_el_ano_de_partida)_5_2001_2006_2011_2014
2
1986	2010
#PARAMETROS_DE_CRECIMIENTO_Y_MORTALIDAD_NATURAL_edades_de_ref_(t1_y_t2)
#_____________________________________________________________________________________
#Loo(17.55)_y_k(1.06)_y_Lo(2.85)_y_alfa(1.357)_y_beta_y_M(1.1)_y_h_(1_linea:parametro_2_linea:CV_3_linea:Fase_Estimacion)
17.55	1.06	0.01	1.357	1e-3	1.1	0.95
0.2	0.2	0.2	0.2	0.2	0.3	0.2
-4	-4	-5	-5	-5	-5	-5
#
#A50_y_rango_y_cv_respectivo
1.0	0.2#Prior_selectividad
0.5	0.2#Sigma_selectividad
#_________________________________________________________________ 
#OPCIONES_Y_FASES_DE_ESTIMACION_SELECTIVIDAD
#(opt_devR)_fase_desviacion_reclutas
1
#(opt_Scru1)_fase_selectividad_crucero 
3
#(opt1_fase)_fase_selectividad_flota 
3
#(opt_q)_fase_de_estimacion_de_q_general
4
#(faseF)_fase_de_estimacion_de_F_general
3
#OTRAS_OPCIONES
#Gama_para_otros_modelos_stock_reclutamiento
1.18
#Opcion_de_Beverton_Holt_(1)_or_Ricker_(2)_or_Desvios_(3)_(opt_Reclu)
3
#CV_prior_F_(la_prior_es_el_vector_de_F_de_la_ultima_evaluacion_*.pin_si_no_hay_prior_relajar_cv_al_maximo)
0.9
#valor_de_PBR_al_RMS
0.86
#Semestres_a_simular_en_el_futuro_(ntime_sim)=nrec_proy_(reclutas_proyectados_fila_84)
4
#Numero_de_corridas(nrun)
1000
#CBA_para_un_semestre(1)_o_para_un_year(2)_(sem_year==)
2
#los_siguientes_2_semestres_a_proyectar_para_la_CBA_(sem1_y_sem2):5_6_opcion_recluta_estable_(3_years_mas)
3	4
#Escenario_de_reclutamiento_(proporcion_de_Ro; Estados de la naturaleza)
1.0
#0.62
#No_bloques_qmph_(primero_es_year_de_partida)
3
1986	1997	2007
#(opt_qmph)_fase_estimacion_q_mph
3
#Sel a la talla para crucero acustico Chile
#en log init_number low_Len50_S;   init_number upp_Len50_S;   init_int phase_Len50_S;
-2	2.5	1
-10.0	5	1
#Opcion_reclutamientos_futuros_estimados
#Opt_reclu_prom_1sem(29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65)
#Opt_reclu_prom_2sem(30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66)
#Opt_reclu_posi_1sem(33,39,47,53,55,57,59,61)
#Opt_reclu_posi_2sem(34,40,48,54,56,58,60,62)
#Opt_reclu_nega_1sem(29,31,35,37,41,43,45,49,51,63,65)
#Opt_reclu_nega_2sem(30,32,36,38,42,44,46,50,52,64,66)
33	39	47	53	55	57	59	61
34	40	48	54	56	58	60	62
