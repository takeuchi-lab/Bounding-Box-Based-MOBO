library(MASS)
library(lpSolve)
library(nsga2R)

seed<-1
total_seed_num<-100
iteration_num<-500
PHV_sample_num<-10000
AF_num<-3			#1:MVA, 2:Proposed, 3:SABBa

simulator_setting<- 1


obj_num<-c(1,2)
obj_dim<-2
input_dim<-6
grid_num<-7

x1_min<- -1
x2_min<- -1
x3_min<- -1
w1_min<- -1
w2_min<- -1
w3_min<- -1

x1_max<- 1
x2_max<- 1
x3_max<- 1
w1_max<- 1
w2_max<- 1
w3_max<- 1

all_x1<- seq(x1_min,x1_max,length=grid_num)
all_x2<- seq(x2_min,x2_max,length=grid_num)
all_x3<- seq(x3_min,x3_max,length=grid_num)
all_w1<- seq(w1_min,w1_max,length=grid_num)
all_w2<- seq(w2_min,w2_max,length=grid_num)
all_w3<- seq(w3_min,w3_max,length=grid_num)

standard_normal_pdf<-function(x)
{
(2*pi)^(-1/2)*exp(-x^2/2)
}

prob_w1<- standard_normal_pdf(all_w1)/sum(standard_normal_pdf(all_w1))
prob_w2<- standard_normal_pdf(all_w2)/sum(standard_normal_pdf(all_w2))
prob_w3<- standard_normal_pdf(all_w3)/sum(standard_normal_pdf(all_w3))


all_x1x2x3<- cbind(cbind(all_x1%x%rep(1,grid_num),rep(1,grid_num)%x%all_x2)%x%rep(1,grid_num),rep(1,grid_num^2)%x%all_x3)
all_w1w2w3<- cbind(cbind(all_w1%x%rep(1,grid_num),rep(1,grid_num)%x%all_w2)%x%rep(1,grid_num),rep(1,grid_num^2)%x%all_w3)
all_w1w2w3_prob<- apply(cbind(cbind(prob_w1%x%rep(1,grid_num),rep(1,grid_num)%x%prob_w2)%x%rep(1,grid_num),rep(1,grid_num^2)%x%prob_w3),1,prod)

all_x1x2x3w1w2w3<-   cbind(all_x1x2x3%x%rep(1,grid_num^3),rep(1,grid_num^3)%x%all_w1w2w3)


sigmaf<-1
sigman<-10^(-6)
length_scale<-4
beta<-3

true_risk_sample<-1000



ParEGO_rho<-0.05
ParEGO_length<-length_scale[1]
ParEGO_sigmaf<-1
ParEGO_sigman<-10^(-8)

SMSego_sample_num<-50

EHVI_sample_num<-20

MESMO_sample_num<-10
RFM_sample_num<-500

PFES_sample_num<-10
max_PF_num<-50


EMMDI_sample_num<-20
JES_sample_num<-10

SABBa_AF_sample_num<-20

SABBa_new_input_sample<-10
SABBa_initial<-2#0.5
SABBa_ratio<-0.9
SABBa_min_ratio<-0.001

###############################################################
##################     True function    #######################
###############################################################

true_function<-function(a)
{
#x1<-a[1]
#x3<-a[2]
#x5<-a[3]
#x2<-a[4]
#x4<-a[5]
#x6<-a[6]

x3<-a[1]
x4<-a[2]
x5<-a[3]
x1<-a[4]
x2<-a[5]
x6<-a[6]

(-1)*(( 100*(x2-x1^2)^2+(1-x1)^2 )+
( 100*(x3-x2^2)^2+(1-x2)^2 )+
( 100*(x4-x3^2)^2+(1-x3)^2 )+
( 100*(x5-x4^2)^2+(1-x4)^2 )+
( 100*(x6-x5^2)^2+(1-x5)^2 )-273.45)/sqrt(28153.22)
}


set.seed(1)
true_risk_sample_w1<-colSums(rmultinom(1000,c(1:grid_num),prob_w1)*all_w1)
true_risk_sample_w2<-colSums(rmultinom(1000,c(1:grid_num),prob_w2)*all_w2)
true_risk_sample_w3<-colSums(rmultinom(1000,c(1:grid_num),prob_w3)*all_w3)



true_risk_function<-function(x)
{
x1<-x[1]
x2<-x[2]
x3<-x[3]


X<- cbind(x1,x2,x3,all_w1w2w3)
result<- apply(X,1,true_function)
mu<-sum(result*all_w1w2w3_prob)
var<-sum((result-mu)^2*all_w1w2w3_prob)
c(mu,-sqrt(var))
}


###############################################################
###############################################################
###############################################################



###############################################################
##################   Input and Output   #######################
###############################################################

set.seed(1)
x1<- runif(100,x1_min,x1_max)
x2<- runif(100,x2_min,x2_max)
x3<- runif(100,x3_min,x3_max)


X<- cbind(x1,x2,x3)

resu<- apply(all_x1x2x3,1,true_risk_function)
#plot(t(resu)[,1],t(resu)[,2])


ref_point<-apply(resu,1,min)-1
max_point<-apply(resu,1,max)+1
vol_refp_maxp<- prod(max_point-ref_point)

threshold_s_1_2<-apply(resu,1,max)-apply(resu,1,min)
threshold_s1<-threshold_s_1_2[1]
threshold_s2<-threshold_s_1_2[2]


###############################################################
###############################################################
###############################################################



###############################################################
##################       Function       #######################
###############################################################




compute_kernel<-function(a)
{
x<-a[1:6]
y<-a[7:12]
d<- sum((x-y)^2)

sigmaf*exp(-d/length_scale)

}

compute_kernel_SABBa<-function(a)
{
x<-a[1:3]
y<-a[4:6]
d<- sum((x-y)^2)


c(sigmaf*exp(-d/length_scale),sigmaf*exp(-d/length_scale))

}





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







compute_PF<-function(A)
{
S<-A
	if(length(S[,1])==1)
	{
	B<-cbind(S,1)
	} else {
	C<-c(1:length(A[,1]))
	B<-c()
		while(length(C)>0)
		{

			

		n<-C[1]
		a<-A[n,]

			if(length(C)==1)
			{
			B<-rbind(B,c(a,n))
			break
			}

		S<-A[C,]
		AAA<-a-t(S)
		AAAA<-apply(AAA,2,min)
		C<-C[AAAA< 0]
		S<-A[C,]
		AAA<- AAA[,colSums(abs(AAA))!=0]

			if(length(AAA)>0)
			{
			aa<-min(apply(cbind(AAA),2,max))
			} else {
			aa<-1
			}


			
			if(aa>0)
			{
			B<-rbind(B,c(a,n))
			}
		}

	}
B
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
#results <- nsga2R(fn=negative_true_risk_function, varNo=2, objDim=2, lowerBounds=c(x1_min,x2_min), upperBounds=c(x1_max,x2_max),
#popSize=50, tourSize=2, generations=50, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)
#plot(results$objectives)






true_PF<- compute_PF_2dim(t(resu))[,1:obj_dim]

#set.seed(1)
true_PF_cell<- create_cells(true_PF,ref_point,max_point,c(1:length(true_PF[,1])))
true_PHV<- sum(apply(rbind(rbind(true_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(true_PF_cell)[,1:(obj_dim)]),1,prod))

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




for(w in 1:total_seed_num)
{
set.seed((seed-1)*total_seed_num+w)

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


new_X_num<- sample(c(1:grid_num^3),1)
new_W_num<- sample(c(1:grid_num^3),1)



new_input<- c(all_x1x2x3[new_X_num,],all_w1w2w3[new_W_num,])
new_Y<- true_function(new_input)+rnorm(1,0,sqrt(sigman))
K<- matrix(numeric(iteration_num^2),iteration_num,iteration_num)	

all_k_vec_mat<- matrix(numeric(grid_num^6*iteration_num),grid_num^6,iteration_num)
	if(AF_num==3)
	{
	SABBa_new_input_counter<-0
	SABBa_new_input_counter_old<-0


	est_mu_tilde_mat<-matrix(numeric(iteration_num*2),iteration_num)
	est_mu_tilde_var_mat<-matrix(numeric(iteration_num*2),iteration_num)
	#kvec_SABBa_1<- matrix(numeric(SABBa_sample_num*iteration_num),SABBa_sample_num,iteration_num)
	#kvec_SABBa_2<- matrix(numeric(SABBa_sample_num*iteration_num),SABBa_sample_num,iteration_num)
	X_tilde<-c(new_input[1:3])
	ratio<-max(SABBa_initial,SABBa_min_ratio)

	}

	for(ww in 1:iteration_num)
	{
	input<-rbind(input,new_input)	
	X_num<-c(X_num,new_X_num)
	W_num<-c(W_num,new_W_num)

	Y<-rbind(Y,new_Y)

	k_vec<-apply(cbind(new_input[1],new_input[2],new_input[3],new_input[4],new_input[5],new_input[6],input),1,compute_kernel)	

	K[ww,1:ww]<- k_vec


		if(ww==1)
		{
		G<- K[1:ww,1:ww]

		} else {
			G<- K[1:ww,1:ww]+t(K[1:ww,1:ww])-diag(diag(K[1:ww,1:ww]))

		}
	C<- G+ sigman*diag(ww)


	C_inv<- solve(C)


	
	C_inv_Y<- as.vector(C_inv%*%Y)

	all_k_vec_mat[,ww]<- apply(cbind(new_input[1],new_input[2],new_input[3],new_input[4],new_input[5],new_input[6],all_x1x2x3w1w2w3),1,compute_kernel)	
	
	post_mean_function<-function(x)
	{
	A<- cbind(x[1],x[2],x[3],x[4],x[5],x[6],input)
	B<- apply(A,1,compute_kernel)
	sum(C_inv_Y*B)

	}

	post_var_function<-function(x)
	{
	A<- cbind(x[1],x[2],x[3],x[4],x[5],x[6],input)
	B<- apply(A,1,compute_kernel)
	as.vector(abs(sigmaf-t(B)%*%C_inv%*%B))
	}

	ucb_function<-function(x)
	{
	post_mean_function(x)+beta*sqrt(post_var_function(x))
	}

	lcb_function<-function(x)
	{
	post_mean_function(x)-beta*sqrt(post_var_function(x))
	}


	
	predicted_PFn<-   compute_PF_2dim(t(apply(rbind(input[,1:3]),1,true_risk_function)))
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
		



		if(AF_num==1)
		{
		starttime<-Sys.time()

		#set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		#ucb_lcb_sample_w1<- runif(100,w1_min,w1_max)
		#ucb_lcb_sample_w2<- runif(100,w2_min,w2_max)



		ucb_lcb_tilde<-function(x)
		{
		A<- cbind(x[1],x[2],x[3],all_w1w2w3)
		ucb_all<- apply(A,1,ucb_function)
		lcb_all<- apply(A,1,lcb_function)

		ucb_F1<- sum(ucb_all*all_w1w2w3_prob)
		lcb_F1<- sum(lcb_all*all_w1w2w3_prob)

		l_tilde<- lcb_all-ucb_F1
		u_tilde<- ucb_all-lcb_F1

		l_tilde2<-l_tilde^2
		u_tilde2<-u_tilde^2

		u_tilde_sq_all<- apply(cbind(l_tilde2,u_tilde2),1,max)
		l_tilde_sq_all<- apply(cbind(l_tilde2,u_tilde2),1,min)*(1-(l_tilde<=0)*(0<=u_tilde))

		ucb_F2<- -sqrt(abs(sum(l_tilde_sq_all*all_w1w2w3_prob)))
		lcb_F2<- -sqrt(abs(sum(u_tilde_sq_all*all_w1w2w3_prob)))

		c(lcb_F1,ucb_F1,lcb_F2,ucb_F2)
		}

		LCB_UCB<- t(apply(all_x1x2x3,1,ucb_lcb_tilde))


		LCB_mat<-LCB_UCB[,c(1,3)]
		UCB_mat<-LCB_UCB[,c(2,4)]


		lcb_PFn<- compute_PF_2dim(LCB_mat)
		lcb_PF<-rbind(lcb_PFn[,1:2])
		estimated_x_num<- lcb_PFn[,3]

		AF_x<-function(i)
		{
		ucb<- UCB_mat[i,]
		lcb<- LCB_mat[i,]

		distX<- min(apply((ucb-t(lcb_PF)),2,max))
		sum((ucb-lcb)^2)*(distX>=0)
		}

		AF_X<- apply( cbind(c(1:grid_num^3)),1,AF_x)
		new_X_num<- order(AF_X,decreasing=T)[1]

		if(simulator_setting==1)
		{
		
		AF_W<-apply(cbind(all_x1x2x3[new_X_num,1]	,all_x1x2x3[new_X_num,2],all_x1x2x3[new_X_num,3],all_w1w2w3),1,post_var_function)
		new_W_num<- order(AF_W,decreasing=T)[1]



		} else {
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		new_W_num<- c(1:grid_num^3)[rmultinom(1,c(1:grid_num^3),all_w1w2w3_prob)==1]
		}


		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		


		inference_predicted_PF_PFn<-compute_PF_2dim(t(apply(rbind(all_x1x2x3[estimated_x_num,]),1,true_risk_function)))
		inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


		inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
		inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


		inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
		inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

		inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
		inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

		inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
		inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]


		}

		
		if(AF_num==2)
		{
		starttime<-Sys.time()

		#set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		#ucb_lcb_sample_w1<- runif(100,w1_min,w1_max)
		#ucb_lcb_sample_w2<- runif(100,w2_min,w2_max)



		ucb_lcb_tilde<-function(x)
		{
		A<- cbind(x[1],x[2],x[3],all_w1w2w3)
		ucb_all<- apply(A,1,ucb_function)
		lcb_all<- apply(A,1,lcb_function)

		ucb_F1<- sum(ucb_all*all_w1w2w3_prob)
		lcb_F1<- sum(lcb_all*all_w1w2w3_prob)

		l_tilde<- lcb_all-ucb_F1
		u_tilde<- ucb_all-lcb_F1

		l_tilde2<-l_tilde^2
		u_tilde2<-u_tilde^2

		u_tilde_sq_all<- apply(cbind(l_tilde2,u_tilde2),1,max)
		l_tilde_sq_all<- apply(cbind(l_tilde2,u_tilde2),1,min)*(1-(l_tilde<=0)*(0<=u_tilde))

		ucb_F2<- -sqrt(abs(sum(l_tilde_sq_all*all_w1w2w3_prob)))
		lcb_F2<- -sqrt(abs(sum(u_tilde_sq_all*all_w1w2w3_prob)))

		c(lcb_F1,ucb_F1,lcb_F2,ucb_F2)
		}

		LCB_UCB<- t(apply(all_x1x2x3,1,ucb_lcb_tilde))


		LCB_mat<-LCB_UCB[,c(1,3)]
		UCB_mat<-LCB_UCB[,c(2,4)]


		lcb_PFn<- compute_PF_2dim(LCB_mat)
		lcb_PF<-rbind(lcb_PFn[,1:2])
		estimated_x_num<- lcb_PFn[,3]

		AF_x<-function(i)
		{
		ucb<- UCB_mat[i,]
		lcb<- LCB_mat[i,]

		 min(apply((ucb-t(lcb_PF)),2,max))
		

		}

		AF_X<- apply( cbind(c(1:grid_num^3)),1,AF_x)
		new_X_num<- order(AF_X,decreasing=T)[1]

		if(simulator_setting==1)
		{
		
		AF_W<-apply(cbind(all_x1x2x3[new_X_num,1]	,all_x1x2x3[new_X_num,2],all_x1x2x3[new_X_num,3],all_w1w2w3),1,post_var_function)
		new_W_num<- order(AF_W,decreasing=T)[1]



		} else {
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		new_W_num<- c(1:grid_num^3)[rmultinom(1,c(1:grid_num^3),all_w1w2w3_prob)==1]
		}


		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		


		inference_predicted_PF_PFn<-compute_PF_2dim(t(apply(rbind(all_x1x2x3[estimated_x_num,]),1,true_risk_function)))
		inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


		inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
		inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


		inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
		inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

		inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
		inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

		inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
		inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]


		}





		if(AF_num==3)
		{
		
		starttime<-Sys.time()
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)

			if(SABBa_new_input_counter==0)
			{

			SABBa_new_input_X<-  all_x1x2x3[sample(c(1:grid_num^3),SABBa_new_input_sample),]#  cbind(runif(SABBa_new_input_sample,x1_min,x1_max), runif(SABBa_new_input_sample,x2_min,x2_max))
			if(ww>1)
			{
			ratio<-  max(SABBa_ratio*ratio,SABBa_min_ratio)
			X_tilde<-rbind(X_tilde[dist_candidate_X>0,])

			}

			X_tilde<-rbind(X_tilde,SABBa_new_input_X)

			BB<-rbind(X_tilde[1,])
			if(length(X_tilde[,1])>1)
			{
				for(abc in 2:length(X_tilde[,1]))
				{
					if(min(colSums(abs(t(BB)-X_tilde[abc,])))>0)
					{
					BB<- rbind(BB,X_tilde[abc,])
					}
				}
			}
			X_tilde<-BB
			X_tilde_current_id<- X_tilde[,1]*0
			for(abc in 1:length(X_tilde[,1]))
			{
			X_tilde_current_id[abc]<-(c(1:grid_num^3)[colSums(abs(t(all_x1x2x3)-X_tilde[abc,]))==0])[1]
			}
			X_tilde_current_id<-sort(X_tilde_current_id)
			prior_SA_var<-  t(apply(cbind(X_tilde,X_tilde),1,compute_kernel_SABBa))
			
			threshold_s1_bar<- threshold_s1*ratio
			threshold_s2_bar<- threshold_s2*ratio
			SABBa_new_input_counter<-1
			}






		



		mu_tilde<-function(x)
		{
		A<- cbind(x[1],x[2],x[3],all_w1w2w3)
		B<- apply(A,1,post_mean_function)
		mu<- sum(B*all_w1w2w3_prob)
		var<-sum((B-mu)^2*all_w1w2w3_prob)
		c(mu,-sqrt(abs(var)))
		}

		error_x_w<-function(x)
		{
		A<- cbind(x[1],x[2],x[3],all_w1w2w3)
		B<- sqrt(apply(A,1,post_var_function))
		epsilon_bar_q_x<- 3*B
		epsilon_bar_mu<-sum(epsilon_bar_q_x*all_w1w2w3_prob)
		q_hat_x<- apply(A,1,post_mean_function)
		mu_tilde<- sum(q_hat_x*all_w1w2w3_prob)
		E1<- epsilon_bar_mu
		E2<- sum(((epsilon_bar_mu+epsilon_bar_q_x)^2+2*abs(q_hat_x-mu_tilde)*(epsilon_bar_mu+epsilon_bar_q_x))*all_w1w2w3_prob)
		c(E1,sqrt(abs(E2)))
		}

	X_tilde_leave_id<-(X_tilde_current_id%x%rep(1,grid_num^3)-1)*grid_num^3+rep(1,length(X_tilde_current_id))%x%c(1:grid_num^3)

	post_mean_mat<- cbind(all_k_vec_mat[X_tilde_leave_id,1:ww])%*%C_inv_Y
	post_var_mat<-   abs(sigmaf -rowSums((cbind(all_k_vec_mat[X_tilde_leave_id,1:ww])%*%C_inv)*cbind(all_k_vec_mat[X_tilde_leave_id,1:ww])))




		mu_tilde_mat<-   colSums(matrix(post_mean_mat,grid_num^3,length(X_tilde_current_id))*all_w1w2w3_prob)
		sd_tilde_mat<-   -sqrt(abs(colSums((matrix(post_mean_mat,grid_num^3,length(X_tilde_current_id))-t(mu_tilde_mat)%x%rep(1,grid_num^3))^2*all_w1w2w3_prob)))
		
		epsilon_bar_q_x<-3*sqrt(abs(matrix(post_var_mat,grid_num^3,length(X_tilde_current_id))))
		epsilon_bar_mu<-colSums(epsilon_bar_q_x*all_w1w2w3_prob)
		q_hat_x<- matrix(post_mean_mat,grid_num^3,length(X_tilde_current_id))
		E1<- epsilon_bar_mu
		E2<- colSums(((t(epsilon_bar_mu)%x%rep(1,grid_num^3)+epsilon_bar_q_x)^2+2*abs(q_hat_x-t(mu_tilde_mat)%x%rep(1,grid_num^3))*(t(epsilon_bar_mu)%x%rep(1,grid_num^3)+epsilon_bar_q_x))*all_w1w2w3_prob)


		Risk_mat<- cbind(mu_tilde_mat,sd_tilde_mat)
		Risk_mat_error<- cbind(E1,sqrt(abs(E2)))

		candidate_X_mu_tilde<- Risk_mat#[X_tilde_current_id,]



		candidate_X_error_x_w<- Risk_mat_error#[X_tilde_current_id,]






		est_mu_tilde_mat[ww,]<-   mu_tilde(input[ww,])
		est_mu_tilde_var_mat[ww,]<-   (error_x_w(input[ww,])/3)^2

		kvec_SABBa<- t(apply(cbind(X_tilde%x%rep(1,ww),rep(1,length(X_tilde[,1]))%x%rbind(input)[1:ww,1],rep(1,length(X_tilde[,1]))%x%rbind(input)[1:ww,2],rep(1,length(X_tilde[,1]))%x%rbind(input)[1:ww,3]),1,compute_kernel_SABBa))

		kvec_SABBa_1<-  t(matrix(kvec_SABBa[,1],ww,length(X_tilde[,1])))
		kvec_SABBa_2<-  t(matrix(kvec_SABBa[,2],ww,length(X_tilde[,1])))

		K_SABBa<-   apply(cbind(rbind(input[,1:3])%x%(numeric(ww)+1),(numeric(ww)+1)%x%rbind(input[,1:3])),1,compute_kernel_SABBa)
		K_SABBa_1<-  matrix(K_SABBa[1,],ww,ww)
		K_SABBa_2<-  matrix(K_SABBa[2,],ww,ww)

		if(ww==1)
		{
		C_SABBa_1<- K_SABBa_1+(est_mu_tilde_var_mat[1:ww,1])
		C_SABBa_2<- K_SABBa_2+(est_mu_tilde_var_mat[1:ww,2])
		} else {
		C_SABBa_1<- K_SABBa_1+diag(est_mu_tilde_var_mat[1:ww,1])
		C_SABBa_2<- K_SABBa_2+diag(est_mu_tilde_var_mat[1:ww,2])
		}
		
		inv_C_SABBa_1<- solve(C_SABBa_1)
		inv_C_SABBa_2<- solve(C_SABBa_2)

		mu_SA_1<-  as.vector(cbind(kvec_SABBa_1[,1:ww])%*%inv_C_SABBa_1%*%est_mu_tilde_mat[1:ww,1])
		mu_SA_2<-  as.vector(cbind(kvec_SABBa_2[,1:ww])%*%inv_C_SABBa_2%*%est_mu_tilde_mat[1:ww,2])

		var_SA_1<- abs(prior_SA_var[,1]-rowSums((cbind(kvec_SABBa_1[,1:ww])%*%inv_C_SABBa_1)*cbind(kvec_SABBa_1[,1:ww])))
		var_SA_2<- abs(prior_SA_var[,2]-rowSums((cbind(kvec_SABBa_2[,1:ww])%*%inv_C_SABBa_2)*cbind(kvec_SABBa_2[,1:ww])))

		error_SA_1<- 3*sqrt(var_SA_1)
		error_SA_2<- 3*sqrt(var_SA_2)

		
		if(SABBa_new_input_counter_old==1)
		{
		X_tilde_num<- length(X_tilde[,1])

			if(sum((1*(error_SA_1[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num]<threshold_s1_bar)*(error_SA_2[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num]<threshold_s2_bar)))>0)
			{
			mu_SA_new<-rbind(mu_SA_old,rbind(cbind(mu_SA_1[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num],mu_SA_2[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num])[(1*(error_SA_1[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num]<threshold_s1_bar)*(error_SA_2[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num]<threshold_s2_bar))==1,]))
			X_tilde_new<- rbind(X_tilde_old, rbind(SABBa_new_input_X[(1*(error_SA_1[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num]<threshold_s1_bar)*(error_SA_2[(X_tilde_num-SABBa_new_input_sample+1):X_tilde_num]<threshold_s2_bar))==1,]))
			

			inference_predicted_PF_PFn<-compute_PF(mu_SA_new)

		reported_Y<-rbind(t(apply(rbind((rbind(X_tilde_new))[inference_predicted_PF_PFn[,(obj_dim+1)],]),1,true_risk_function)))
		
		



		inference_predicted_PF_PFn<-compute_PF_2dim(rbind(t(apply(rbind((rbind(X_tilde_new))[inference_predicted_PF_PFn[,(obj_dim+1)],]),1,true_risk_function))))
		inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


		inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
		inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


		#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
		inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

		#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
		inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

		#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
		inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y),0)



			}

		SABBa_new_input_counter_old<-0
		}





		lcb_pf<-rbind(compute_PF_2dim(candidate_X_mu_tilde-candidate_X_error_x_w)[,1:obj_dim])
		dist_candidate_X<-compute_dist(candidate_X_mu_tilde+candidate_X_error_x_w,lcb_pf)

		check<-(dist_candidate_X>0)*(1-(error_SA_1<threshold_s1_bar)*(error_SA_2<threshold_s2_bar))*(1-(candidate_X_error_x_w[,1]<threshold_s1_bar)*(candidate_X_error_x_w[,2]<threshold_s2_bar))
				

		if(sum(check)==0)
		{
		new_input_x<- X_tilde[sample(c(1:length(X_tilde[,1])),1),]
		new_X_num<-c(1:grid_num^3)[colSums(abs(t(all_x1x2x3)-new_input_x))==0]

		SABBa_new_input_counter<-0

		} else {
			if(sum(check)==1)
			{
			new_input_x<- X_tilde[check==1,]
			new_X_num<-c(1:grid_num^3)[colSums(abs(t(all_x1x2x3)-new_input_x))==0]

			} else {
			B_minus<- candidate_X_mu_tilde-candidate_X_error_x_w
			B_plus<- candidate_X_mu_tilde+candidate_X_error_x_w




				function_pop<-function(i)
				{
				n<- length(B_minus[,1])-1
				v<- setdiff(c(1:length(B_minus[,1])),i)
				A<- cbind(runif(SABBa_AF_sample_num,B_minus[i,1],B_plus[i,1]),runif(SABBa_AF_sample_num,B_minus[i,2],B_plus[i,2])	)
				result<- numeric(n)
					for(j in 1:n)
					{
					B<- cbind(runif(SABBa_AF_sample_num,B_minus[v[j],1],B_plus[v[j],1]),runif(SABBa_AF_sample_num,B_minus[v[j],2],B_plus[v[j],2])	)
					result[j]<-mean(1-apply(B-A,1,min)>0)
					}
				min(result)
				}
			AF<- apply( cbind(c(1:length(X_tilde[,1]))),1,function_pop)
			AF<-AF*check
			new_input_x<-X_tilde[order(AF,decreasing=T)[1],]
			new_X_num<-c(1:grid_num^3)[colSums(abs(t(all_x1x2x3)-new_input_x))==0]

			}
		}

		if(simulator_setting==1)
		{
			

		
		




			A<- cbind(new_input_x[1],new_input_x[2],new_input_x[3],all_w1w2w3)
			B<- sqrt(apply(A,1,post_var_function))
			epsilon_bar_q_x<- 3*B
			ck_xi<-epsilon_bar_q_x*all_w1w2w3_prob
			

			c_bar_k<- (ck_xi-min(ck_xi))/(max(ck_xi)-min(ck_xi))

			new_X_tilde_num<-			c(1:length(rbind(X_tilde)[,1]))[colSums(abs(t(rbind(X_tilde))-new_input_x))==0]
			new_X_tilde_num<-new_X_tilde_num[1]
			wei<-(candidate_X_error_x_w[new_X_tilde_num,]/threshold_s2_bar)^2
			
		
			
	              AF_w<- rowSums(c_bar_k%x%t(wei))
			new_W_num<- order(AF_w,decreasing=T)[1]



		} else {
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		new_W_num<- c(1:grid_num^3)[rmultinom(1,c(1:grid_num^3),all_w1w2w3_prob)==1]
		}


		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 



		if(SABBa_new_input_counter==0)
		{
		
		inference_predicted_PF_PFn<-compute_PF_2dim(cbind(mu_SA_1[check==0],mu_SA_2[check==0]))

		reported_Y_mat<-rbind(t(apply(rbind((rbind(X_tilde[check==0,]))[inference_predicted_PF_PFn[,(obj_dim+1)],]),1,true_risk_function)))

		mu_SA_old<-cbind(mu_SA_1[check==0],mu_SA_2[check==0])
		X_tilde_old<- rbind(X_tilde[check==0,])
		
		SABBa_new_input_counter_old<-1



		inference_predicted_PF_PFn<-compute_PF_2dim(rbind(t(apply(rbind((rbind(X_tilde[check==0,]))[inference_predicted_PF_PFn[,(obj_dim+1)],]),1,true_risk_function))))
		inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


		inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
		inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


		#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
		inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

		#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
		inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

		#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
		inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y_mat),0)
		} else {
			if(ww==1)
			{
			inference_predicted_PF_PFn<-compute_PF_2dim(cbind(mu_SA_1,mu_SA_2))

			reported_Y_mat<-rbind(t(apply(rbind((rbind(X_tilde))[inference_predicted_PF_PFn[,(obj_dim+1)],]),1,true_risk_function)))

			inference_predicted_PF_PFn<-compute_PF_2dim(rbind(t(apply(rbind((rbind(X_tilde))[inference_predicted_PF_PFn[,(obj_dim+1)],]),1,true_risk_function))))
			inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


			inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
			inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


			#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
			inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

			#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
			inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

			#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
			inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y_mat),0)
			} else {
			inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret_mat[w,(ww-1)]

			inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy_mat[w,(ww-1)]

			inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2_mat[w,(ww-1)]


			}

			}


		

		}



	######################################################################
	set.seed(((seed-1)*total_seed_num+w)*1000+ww)
	new_input<- c(all_x1x2x3[new_X_num,],all_w1w2w3[new_W_num,])
	new_Y<- true_function(new_input)+rnorm(1,0,sqrt(sigman))
		
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




