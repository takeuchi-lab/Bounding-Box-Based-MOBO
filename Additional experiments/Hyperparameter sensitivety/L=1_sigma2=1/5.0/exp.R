library(MASS)
library(lpSolve)
library(nsga2R)

length_scale<-1
sigmaf<-1


seed<-1
total_seed_num<-100
iteration_num<-500

epsilon_PF<-0.0
epsilon_ucb<-0.0
epsilon_lcb<-0.0
epsilon_AF_x<- 0.0
epsilon_AF_w<- 0.0

main_name<- "GP sample path"

simulator_setting<- 1


obj_num<-c(1,2)
obj_dim<-2
input_dim<-4
grid_num<-16

x1_min<- -2
x2_min<- -2
w1_min<- -2



x1_max<- 2
x2_max<- 2
w1_max<- 2



all_x1<- seq(x1_min,x1_max,length=grid_num)
all_x2<- seq(x2_min,x2_max,length=grid_num)

all_w1<- seq(w1_min,w1_max,length=grid_num)



standard_normal_pdf<-function(x)
{
(2*pi)^(-1/2)*exp(-x^2/2)
}

prob_w1<- standard_normal_pdf(all_w1)/sum(standard_normal_pdf(all_w1))




all_x1x2<- cbind(all_x1%x%rep(1,grid_num),rep(1,grid_num)%x%all_x2)

all_x1x2w1<- cbind(all_x1x2%x%rep(1,grid_num),rep(1,grid_num^2)%x%all_w1)
all_w1_prob<- prob_w1



sigman<-10^(-6)





###############################################################
##################     True function    #######################
###############################################################



compute_kernel<-function(a)
{
x<-a[1:3]
y<-a[4:6]
d<- sum((x-y)^2)

sigmaf*exp(-d/length_scale)

}


compute_kernel2<-function(a1,a2,a3,a4,a5,a6)
{


d<- (a1-a4)^2+(a2-a5)^2+(a3-a6)^2

sigmaf*exp(-d/length_scale)

}



gram_mat<-diag(grid_num^3)
for(i in 1:grid_num^3)
{
gram_mat[i,]<- compute_kernel2(all_x1x2w1[i,1],all_x1x2w1[i,2],all_x1x2w1[i,3],all_x1x2w1[,1],all_x1x2w1[,2],all_x1x2w1[,3])
}

###############################################################
###############################################################
###############################################################


###############################################################
##################   Input and Output   #######################
###############################################################











compute_PF_2dim<-function(A)
{
S<-A
l<- c(1:length(S[,1]))	
n<-c()
C<-c()

	while(length(l)>0)	
	{
	m2<-max(S[,2])
	m1<-max(S[S[,2]==m2,1])
	C<-rbind(C,c(m1,m2))
	
	n<-c(n,l[(S[,1]==m1)&(S[,2]==m2)][1])
	l<-l[S[,1]>m1]	
	S<-rbind(A[l,])
	}
B<-cbind(C,n)
rbind(B[order(B[,3]),])
}






compute_PHV<-function(A)
{
U<-matrix(numeric(PHV_sample_num*obj_dim),PHV_sample_num,obj_dim)
	for(i in 1:obj_dim)
	{
	U[,i]<- runif(PHV_sample_num,ref_point[i],max_point[i])
	}

a<-0
	for(i in 1:PHV_sample_num)
	{
	a<-a+(max(apply((t(A)-U[i,]),2,min))>0)*1
	}
(a/PHV_sample_num)*vol_refp_maxp
}


compute_e_accuracy<-function(truePF,predPF)
{
d<-length(truePF[,1])
a<-numeric(d)
	for(i in 1:d)
	{
	a[i]<-min(apply((truePF[i,]-t(predPF)),2,max))
	}
max(a)
}

compute_e_accuracy2<-function(truePF,predPF)
{
d<-length(truePF[,1])
a<-numeric(d)
	for(i in 1:d)
	{
	a[i]<-max(apply((truePF[i,]-t(predPF)),2,min))
	}
max(a)
}


compute_dist<-function(UCB,LCB)
{
d<-length(UCB[,1])
a<-numeric(d)
	for(i in 1:d)
	{
	a[i]<-min(apply((UCB[i,]-t(LCB)),2,max))
	}
a
}

compute_SMSego_PHVI<-function(y,PF,n)
{
PF1<- rbind(PF)
A<-t(PF1)-y
	if(max(apply(A,2,min))>=0)
	{
	b<-0
	} else {
	B<-t(matrix(runif(n*obj_dim,apply(cbind(ref_point,y),1,min),apply(cbind(max_point,y),1,max)),obj_dim,n))
	B1<- rbind(B[apply((y- t(B)),2,min)>=0,])
	b<-numeric(length(B1[,1]))
		if(length(b)==0)
		{
		b<-0
		} else {
			for(i in 1:length(B1[,1]))
			{
			b[i]<-  (max(apply((t(PF1)-B1[i,]),2,min))<0)*1
			}
			b<-sum(b)/n*prod(apply(cbind(max_point,y),1,max)-apply(cbind(ref_point,y),1,min))
		}
	}
b
}


compute_MESMO<-function(PF,post_mean,post_var)
{

y_j_s_star<-apply(PF,2,max)
gamma_j_s_star<- (y_j_s_star-post_mean)/sqrt(post_var)
result<-sum((gamma_j_s_star*standard_normal_pdf(gamma_j_s_star))/(2*pnorm(gamma_j_s_star))-log(pnorm(gamma_j_s_star)))
result[is.na(result)]<-0
result
}








dummy_mat<-c(0,1)
for(i in 2:obj_dim)
{
ones_v<-numeric(2^(i-1))+1
dummy_mat<- cbind(dummy_mat%x%c(1,1),ones_v%x%c(0,1))

}
dummy_mat<-dummy_mat[-c(1,2^(obj_dim)),]

create_cells<-function(mo_set, ref_point, max_point,mo_index)
{
mo_set<-rbind(mo_set)

  N<-length(mo_set[,1])
  L<-length(mo_set[1,])


    
    if(N == 1)
	{
        return(rbind(c(ref_point,as.vector(mo_set))))
	} else{
        hv <- apply((t(mo_set)-ref_point),2,prod)
        pivot_index <- order(hv,decreasing=T)[1]
        pivot <- mo_set[pivot_index,]

        lower <-ref_point
        upper <- pivot

        for (i in 1:(2^L-2))
	{
            iter_index <- dummy_mat[i,]==0

            new_ref <- ref_point
            new_ref[iter_index] <- pivot[iter_index]

            new_ref_max <- pivot
            new_ref_max[iter_index] <- max_point[iter_index]

            dominate_index <- apply((t(mo_set)- new_ref),2,min)>0
            new_mo_index <- mo_index[dominate_index]
            new_mo_set <- rbind(mo_set[dominate_index, ])
		pro_mat<-(t(t(new_mo_set)-new_ref_max)>0)*1
            new_mo_set<- new_mo_set*(1-pro_mat)+pro_mat*((new_mo_set[,1]*0+1)%x%t(new_ref_max))
	            if (length(new_mo_set) > 0)
			{
                child<-create_cells(new_mo_set, new_ref, new_ref_max, new_mo_index)
                lower <- rbind(lower, child[,1:L])
                upper <- rbind(upper, child[,(L+1):(2*L)])

			}
      }
	}
    return( cbind(rbind(lower), rbind(upper)) )
}









compute_PFES<-function(pm,pv,cell)
{
M<-length(cell[,1])
l<-length(cell[1,])/2
Z<-0
ZmGamma<-0




A<-t(c(pm,pm))%x%(numeric(M)+1)
alpha<- t(t((cell-A))/sqrt(pv))
UP<- rbind(alpha[,(l+1):(2*l)])
LO<-rbind(alpha[,1:l])

Zml<-(pnorm(UP)-pnorm(LO))#+10^(-16)
Zm<- apply(Zml,1,prod)
Z<-sum(Zm)
Gammaml<- (LO*standard_normal_pdf(LO)-UP*standard_normal_pdf(UP))/(Zml)
ZmGamma<- sum(Zm*rowSums(Gammaml))/2



  result<- -log(Z) - ZmGamma/Z
result[is.na(result)]<-0
result
}






compute_EMMDI<-function(sample,cPF)
{
d<-length(sample[,1])
a<-numeric(d)
	for(i in 1:d)
	{
	a[i]<-max(min(apply((sample[i,]-t(cPF)),2,max)),0)
	}
mean(a)
}


compute_JES<-function(pv,after_pm,after_pv,cell)
{
first_term<-(obj_dim/2)*log(2*pi*exp(1))+(1/2)*sum(log(pv+sigman))


M<-length(cell[,1])
l<-length(cell[1,])/2
Z<-0
ZmGamma<-0




A<-t(c(after_pm,after_pm))%x%(numeric(M)+1)
alpha<- t(t((cell-A))/sqrt(after_pv))
UP<- rbind(alpha[,(l+1):(2*l)])
LO<-rbind(alpha[,1:l])

Wjm<-(pnorm(UP)-pnorm(LO))
Wj<- apply(Wjm,1,prod)
W<-sum(Wj)
Gjm<-standard_normal_pdf(UP)-standard_normal_pdf(LO)
Vjm<- UP*standard_normal_pdf(UP)-LO*standard_normal_pdf(LO)

Wj_Gjm_invWjm<-Wj*Gjm/Wjm
Gjm_invWjm<-Gjm/Wjm
invWsumWGW<- colSums(Wj_Gjm_invWjm)/W
AAA<- colSums(Wj*Vjm/Wjm)+W*invWsumWGW^2


covmat<- diag(obj_dim)*0
	for(i in 1:(obj_dim-1))
	{
		for(j in (i+1):obj_dim)
		{
		sqrt(after_pv[i]*after_pv[j])*(1/W)*(sum(Wj_Gjm_invWjm[,i]*Gjm_invWjm[,j])-invWsumWGW[j]*sum(Wj_Gjm_invWjm[,i]))
		}
	}

diag_element<- after_pv+sigman-after_pv*AAA/W

covmat2<- covmat+ t(covmat)+diag(diag_element)

second_term<- obj_dim/2*log(2*i*exp(1))*(1/2)*log(det(covmat2))

  result<- first_term-second_term
result[is.na(result)]<-0
result[result==Inf]<-0
result
}


###############################################################
###############################################################
###############################################################


###############################################################
##################     Pareto Front     #######################
###############################################################
#set.seed(1)
#results <- nsga2R(fn=negative_true_Bayes_risk_function, varNo=2, objDim=2, lowerBounds=c(x1_min,x2_min), upperBounds=c(x1_max,x2_max),
#popSize=50, tourSize=2, generations=50, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)
#plot(results$objectives)



###############################################################
###############################################################
###############################################################


###############################################################
##################     Random feature       ###################
###############################################################
#set.seed(1)
#RFM_mat<-  matrix(numeric(total_input_num*RFM_sample_num),total_input_num,RFM_sample_num)%x%numeric(obj_dim)
#for(i in 1:obj_dim)
#{
#RFM_w_mat<-  mvrnorm(RFM_sample_num,numeric(input_dim),(1/(2/length_scale[obj_num[i]]))*diag(input_dim))
#RFM_b_vec<- runif(RFM_sample_num,0,2*pi)
#	for(j in 1:total_input_num)
#	{
#	RFM_mat[((i-1)*total_input_num+j),]<-     sqrt(2*sigmaf[obj_num[i]]/RFM_sample_num)*cos(as.vector(RFM_w_mat%*%all_input[j,])+RFM_b_vec)
#	}
#}



###############################################################
###############################################################
###############################################################


###############################################################
##################      Result mat      #######################
###############################################################

simple_PHV_regret_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
simple_e_accuracy_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
simple_e_accuracy2_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)

inference_PHV_regret_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
inference_e_accuracy_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
inference_e_accuracy2_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)

#X_num_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
AF_computation_time_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)

###############################################################
###############################################################
###############################################################


###############################################################
##################     Optimization     #######################
###############################################################

#Kvec_mat<-matrix(numeric(total_input_num*iteration_num),total_input_num,iteration_num)%x%numeric(obj_dim)
#post_mean_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
#post_var_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)

#after_post_mean_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
#after_post_var_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)

#ucb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
#lcb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
#post_cov_mat<-diag(total_input_num)%x%numeric(obj_dim)

#svd_mat<-diag(total_input_num)%x%numeric(obj_dim)


threshold_1_vec<- numeric(total_seed_num)
threshold_2_vec<- numeric(total_seed_num)

set.seed(1)
#sample_path<-mvrnorm(2,numeric(grid_num^3),gram_mat)

true_function_mat<- cbind(scan("sample_path_f1.txt"),scan("sample_path_f2.txt"))



true_risk_function_mat<-cbind(colSums(matrix(true_function_mat[,1],grid_num,grid_num^2)*all_w1_prob),colSums(matrix(true_function_mat[,2],grid_num,grid_num^2)*all_w1_prob))

pdf("true_risk_functions.pdf")
plot(true_risk_function_mat[,1],true_risk_function_mat[,2],xlab="Risk1",ylab="Risk2",main=main_name,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
dev.off()

for(w in 1:total_seed_num)
{
set.seed((seed-1)*total_seed_num+w)





ref_point<-apply(true_risk_function_mat,2,min)
max_point<-apply(true_risk_function_mat,2,max)
vol_refp_maxp<- prod(max_point-ref_point)








true_PF<- compute_PF_2dim(true_risk_function_mat)[,1:obj_dim]

threshold_1_vec[w]<-epsilon_PF+epsilon_ucb+epsilon_lcb+epsilon_AF_x+epsilon_AF_w




#set.seed(1)
true_PF_cell<- create_cells(true_PF,ref_point,max_point,c(1:length(true_PF[,1])))
true_PHV<- sum(apply(rbind(rbind(true_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(true_PF_cell)[,1:(obj_dim)]),1,prod))


ep_ac_true_PF<- true_PF-(epsilon_PF+epsilon_ucb+epsilon_lcb+epsilon_AF_x+epsilon_AF_w)

ep_ac_true_PF<-cbind(apply(cbind(ep_ac_true_PF[,1],ref_point[1]),1,max),apply(cbind(ep_ac_true_PF[,2],ref_point[2]),1,max))


ep_ac_true_PF_cell<- create_cells(ep_ac_true_PF,ref_point,max_point,c(1:length(ep_ac_true_PF[,1])))
ep_ac_true_PHV<- sum(apply(rbind(rbind(ep_ac_true_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(ep_ac_true_PF_cell)[,1:(obj_dim)]),1,prod))



threshold_2_vec[w]<- true_PHV-ep_ac_true_PHV

input<-c()
Y<-c()

X_num<-c()
W_num<-c()


simple_PHV_regret<-c()
simple_e_accuracy<-c()
simple_e_accuracy2<-c()


inference_PHV_regret<-c()
inference_e_accuracy<-c()
inference_e_accuracy2<-c()


new_X_num<- sample(c(1:grid_num^2),1)
new_W_num<- sample(c(1:grid_num),1)



new_input<- c(all_x1x2[new_X_num,],all_w1[new_W_num])


new_Y<- true_function_mat[((new_X_num-1)*grid_num+new_W_num),]+rnorm(2,0,sqrt(sigman))
K<- matrix(numeric(iteration_num^2),iteration_num,iteration_num)	
k_vec_mat<- numeric(grid_num^3)%*%t(numeric(iteration_num))


	for(ww in 1:iteration_num)
	{
	beta<- 5#sqrt(2*log(2*grid_num^3*pi^2*ww^2/(6*0.1)))

	input<-rbind(input,new_input)	
	X_num<-c(X_num,new_X_num)
	W_num<-c(W_num,new_W_num)

	Y<-rbind(Y,new_Y)

	k_vec<-apply(cbind(new_input[1],new_input[2],new_input[3],input),1,compute_kernel)	

	K[ww,1:ww]<- k_vec

	
	if(ww<=3)
	{
		if(ww==1)
		{
		G<- K[1:ww,1:ww]

		} else {
			G<- K[1:ww,1:ww]+t(K[1:ww,1:ww])-diag(diag(K[1:ww,1:ww]))

		}
	}

	k_vec_mat[,ww]<- 	apply(cbind(all_x1x2w1,new_input[1],new_input[2],new_input[3]),1,compute_kernel)

	if(ww>3)
	{
		current_num<-(new_X_num-1)*grid_num+new_W_num

	new_C_n_plus_one_Y<- (-post_mean_mat[current_num,]+new_Y)/(post_var_mat[current_num,]+sigman)

	C_inv_A_12<- C_inv%*%k_vec[1:(ww-1)]
	
	c_var_noise<-post_var_mat[current_num,1]+sigman
	update_C_inv_11<- C_inv+C_inv_A_12%x%t(C_inv_A_12)/c_var_noise
	update_C_inv_12<- -C_inv_A_12/c_var_noise
	update_C_inv_21<- t(update_C_inv_12)
	update_C_inv_22<- 1/c_var_noise
	updata_C_inv_21_22<-cbind(update_C_inv_21,update_C_inv_22)
	update_C_inv<- rbind(cbind(update_C_inv_11,update_C_inv_12),updata_C_inv_21_22)

	k_vec_mat_1_ww_minus<-	k_vec_mat[,1:(ww-1)]
	k_vec_mat_ww<-k_vec_mat[,ww]
	k_vec_mat_1_ww<- cbind(k_vec_mat_1_ww_minus,k_vec_mat_ww)

	J1<-k_vec_mat_1_ww_minus%*%C_inv_A_12
	S1<-k_vec_mat_ww-J1
	
	pm2<-post_mean_mat+S1%x%t(new_C_n_plus_one_Y)

	update_pv2<-rowSums((J1%x%t(C_inv_A_12)/c_var_noise)*k_vec_mat_1_ww_minus)
	update_pv3<-rowSums((k_vec_mat_1_ww_minus%*%update_C_inv_12)*k_vec_mat_ww)
	update_pv4<-   rowSums((k_vec_mat_ww%x%updata_C_inv_21_22)*k_vec_mat_1_ww)

	update_pv<-post_var_mat[,1]-update_pv2-update_pv3-update_pv4
	
	C_inv<- update_C_inv
	post_mean_mat<- pm2
	post_var_mat<-  cbind(abs(update_pv),abs(update_pv))

	} else {




	C<- G+ sigman*diag(ww)


	C_inv<- solve(C)


	
	C_inv_Y<- C_inv%*%Y



	

	k_vec_C_inv<-cbind(k_vec_mat[,1:ww])%*%C_inv


	post_mean_mat<- k_vec_C_inv%*%Y
	post_var_mat<- abs(sigmaf-rowSums(k_vec_C_inv*cbind(k_vec_mat[,1:ww])))%x%t(c(1,1))
	
	}


	u_mat<- post_mean_mat+beta*sqrt(post_var_mat)
	l_mat<- post_mean_mat-beta*sqrt(post_var_mat)

	selected_Y<- rbind(true_risk_function_mat[X_num,]	)


	predicted_PFn<-   compute_PF_2dim(selected_Y)
	predicted_PF<-rbind(predicted_PFn[,1:2])

	predicted_PF_cell<-create_cells(predicted_PF,ref_point,apply(predicted_PF,2,max),c(1:length(predicted_PF[,1])))

	predicted_PF_PHV<-  sum(apply(rbind(predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-predicted_PF_cell[,(1):(obj_dim)]),1,prod))
	
	
	simple_PHV_regret<- c(simple_PHV_regret,max(true_PHV-predicted_PF_PHV,0))
	simple_PHV_regret_mat[w,ww]<-	min(simple_PHV_regret)

	simple_e_accuracy<- c(simple_e_accuracy,max(compute_e_accuracy(true_PF,predicted_PF),0))
	simple_e_accuracy_mat[w,ww]<-	min(simple_e_accuracy)

	simple_e_accuracy2<- c(simple_e_accuracy2,max(compute_e_accuracy2(true_PF,predicted_PF),0))
	simple_e_accuracy2_mat[w,ww]<-	min(simple_e_accuracy2)


	
	

	
	################## Acquisition function ##############################
		

		


		starttime<-Sys.time()

		#set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		#ucb_lcb_sample_w1<- runif(100,w1_min,w1_max)
		#ucb_lcb_sample_w2<- runif(100,w2_min,w2_max)



		#ucb_lcb_tilde<-function(x)
		#{
		#A<- cbind(x[1],x[2],all_w1w2)
		#ucb<- colSums(t(apply(A,1,ucb_function))*all_w1w2_prob)
		#lcb<- colSums(t(apply(A,1,lcb_function))*all_w1w2_prob)
		#c(lcb[1],ucb[1],lcb[2],ucb[2])
		#}

		#LCB_UCB<- t(apply(all_x1x2,1,ucb_lcb_tilde))
		LCB_mat1<-colSums(matrix(l_mat[,1],grid_num,grid_num^2)*all_w1_prob)
		LCB_mat2<-colSums(matrix(l_mat[,2],grid_num,grid_num^2)*all_w1_prob)
		UCB_mat1<-colSums(matrix(u_mat[,1],grid_num,grid_num^2)*all_w1_prob)
		UCB_mat2<-colSums(matrix(u_mat[,2],grid_num,grid_num^2)*all_w1_prob)


		LCB_mat<-cbind(LCB_mat1,LCB_mat2)- cbind(runif(grid_num^2,-epsilon_lcb,epsilon_lcb),runif(grid_num^2,-epsilon_lcb,epsilon_lcb))
		UCB_mat<-cbind(UCB_mat1,UCB_mat2)+ cbind(runif(grid_num^2,-epsilon_ucb,epsilon_ucb),runif(grid_num^2,-epsilon_ucb,epsilon_ucb))

		UCB_mat[,1]<-apply(cbind(UCB_mat[,1],LCB_mat[,1]),1,max)
		UCB_mat[,2]<-apply(cbind(UCB_mat[,2],LCB_mat[,2]),1,max)


		lcb_PFn<- compute_PF_2dim(LCB_mat)
		lcb_PF<-rbind(lcb_PFn[,1:2])

		dist_lcb_lcbPF<-compute_dist(LCB_mat,lcb_PF)
		dist_lcb_lcbPF_num<-c(1:grid_num^2)[dist_lcb_lcbPF>=(-epsilon_PF)	]
		sort_dist_lcb_lcbPF_num<- dist_lcb_lcbPF_num[order(dist_lcb_lcbPF[dist_lcb_lcbPF_num])]
		lcb_PF2<-rbind(LCB_mat[sort_dist_lcb_lcbPF_num[1],])
		
		lcb_PF2_counter<-1
		while(compute_e_accuracy(lcb_PF,lcb_PF2)>epsilon_PF)
		{
		lcb_PF2_counter<- lcb_PF2_counter+1
		lcb_PF2<- rbind(lcb_PF2,LCB_mat[sort_dist_lcb_lcbPF_num[lcb_PF2_counter],])
		}

		lcb_PF<-lcb_PF2
		estimated_x_num<- sort_dist_lcb_lcbPF_num[1:length(lcb_PF[,1])]


		

		AF_X<- compute_dist(UCB_mat,LCB_mat)
		AF_X_max<- max(AF_X)
		epsilon_AF_X_num<-   c(1:grid_num^2)[AF_X>=AF_X_max-epsilon_AF_x]
		new_X_num<- epsilon_AF_X_num[order(AF_X[epsilon_AF_X_num])[1]]

		if(simulator_setting==1)
		{
		

		AF_W<-rowSums(beta*sqrt(post_var_mat[((new_X_num-1)*grid_num+c(1:grid_num)),]))

		AF_W_max<- max(AF_W)
		epsilon_AF_W_num<-   c(1:grid_num)[AF_W>=AF_W_max-epsilon_AF_w]
		new_W_num<- epsilon_AF_W_num[order(AF_W[epsilon_AF_W_num])[1]]




		} else {
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		new_W_num<- c(1:grid_num)[rmultinom(1,c(1:grid_num),all_w1_prob)==1]
		}


		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		reported_Y<-rbind(true_risk_function_mat[estimated_x_num,])


		inference_predicted_PF_PFn<-compute_PF_2dim(reported_Y)
		inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


		inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
		inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


		inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
		inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

		inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
		inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

		inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,reported_Y),0))
		inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]


		





	######################################################################
	set.seed(((seed-1)*total_seed_num+w)*1000+ww)
	new_input<- c(all_x1x2[new_X_num,],all_w1[new_W_num])


	new_Y<- true_function_mat[((new_X_num-1)*grid_num+new_W_num),]+rnorm(2,0,sqrt(sigman))
		
	}

}


write.table(simple_PHV_regret_mat,"simple_PHV_regret_mat.txt")
write.table(simple_e_accuracy_mat,"simple_e_accuracy_mat.txt")
write.table(simple_e_accuracy2_mat,"simple_e_accuracy2_mat.txt")


write.table(inference_PHV_regret_mat,"inference_PHV_regret_mat.txt")
write.table(inference_e_accuracy_mat,"inference_e_accuracy_mat.txt")
write.table(inference_e_accuracy2_mat,"inference_e_accuracy2_mat.txt")


write.table(AF_computation_time_mat,"AF_computation_time_mat.txt")
#write.table(X_num_mat,"X_num_mat.txt")

write(threshold_1_vec,"threshold_1_vec.txt")
write(threshold_2_vec,"threshold_2_vec.txt")








