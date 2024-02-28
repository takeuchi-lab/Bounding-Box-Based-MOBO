library(MASS)
library(lpSolve)
library(nsga2R)



seed<-1
total_seed_num<-100

naive_sample_num<-5
iteration_num<-150/naive_sample_num

PHV_sample_num<-10000
AF_num<-6			#1:random, 2:US, 3:EHVI, 4:ePAL, 5:proposed, 6:EMmI



obj_num<-c(1,2)
obj_dim<-2
input_dim<-4
x1_min<- 0.25
x2_min<- 0.25
w1_min<- -0.25
w2_min<- -0.25

x1_max<- 0.75
x2_max<- 0.75
w1_max<- 0.25
w2_max<- 0.25


sigmaf<-c(1,1)
sigman<-c(10^(-6),10^(-6))
length_scale<-c(0.2,0.2)
beta<-3

true_Bayes_risk_sample<-1000
naive_sample<-10
nsga2_pop_num<-50

ParEGO_rho<-0.05
ParEGO_length<-length_scale[1]
ParEGO_sigmaf<-2
ParEGO_sigman<-10^(-8)

SMSego_sample_num<-50

EHVI_sample_num<-20

MESMO_sample_num<-10
RFM_sample_num<-500

PFES_sample_num<-10
max_PF_num<-50


EMMDI_sample_num<-20
JES_sample_num<-10

SABBa_sample_num<-1000

###############################################################
##################     True function    #######################
###############################################################

true_function<-function(x)
{
x1<-x[1]
x2<-x[2]
w1<-x[3]
w2<-x[4]

a<-x1+w1
b<-x2+w2

f1<-a
g<- 1+9*b
#h<-1-sqrt(f1/g)
f2<-g-sqrt(abs(f1*g))
#c(-(f1-0.5)/sqrt(0.041),-(f2-24)/sqrt(380))
-c((f1-0.5)/sqrt(0.042),(f2-3.9085)/sqrt(2.5615))
}

true_Bayes_risk_function<-function(x)
{
x1<-x[1]
x2<-x[2]

f1<-x1
f2<- 1+9*x2-(w1_max-w1_min)^(-1)*(2/3)*((abs(w1_max+x1))^(3/2)-(abs(w1_min+x1))^(3/2))*
(w2_max-w2_min)^(-1)*(2/3)*((abs(1+9*(w2_max+x2)))^(3/2)-(abs(1+9*(w2_min+x2)))^(3/2))/9
-c((f1-0.5)/sqrt(0.042),(f2-3.9085)/sqrt(2.5615))
}



negative_true_Bayes_risk_function<-function(x)
{
x1<-x[1]
x2<-x[2]

f1<-x1
f2<- 1+9*x2-(w1_max-w1_min)^(-1)*(2/3)*((abs(w1_max+x1))^(3/2)-(abs(w1_min+x1))^(3/2))*
(w2_max-w2_min)^(-1)*(2/3)*((abs(1+9*(w2_max+x2)))^(3/2)-(abs(1+9*(w2_min+x2)))^(3/2))/9
c((f1-0.5)/sqrt(0.042),(f2-3.9085)/sqrt(2.5615))
}

###############################################################
###############################################################
###############################################################



###############################################################
##################   Input and Output   #######################
###############################################################

set.seed(1)
x1<- runif(1000,x1_min,x1_max)
x2<- runif(1000,x2_min,x2_max)

X<- cbind(x1,x2)

resu<- apply(X,1,true_Bayes_risk_function)



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
standard_normal_pdf<-function(x)
{
(2*pi)^(-1/2)*exp(-x^2/2)
}



compute_kernel<-function(a)
{
x<-a[1:4]
y<-a[5:8]
d<- sum((x-y)^2)

c(sigmaf[1]*exp(-d/length_scale[1]),sigmaf[2]*exp(-d/length_scale[2]))

}


compute_kernel_naive<-function(a)
{
x<-a[1:2]
y<-a[3:4]
d<- sum((x-y)^2)

c(sigmaf[1]*exp(-d/(length_scale[1]/2)),sigmaf[2]*exp(-d/(length_scale[2]/2)))

}



compute_kernel_SABBa<-function(a)
{
x<-a[1:2]
y<-a[3:4]
d<- sum((x-y)^2)

c(sigmaf[1]*exp(-d/length_scale[1]),sigmaf[2]*exp(-d/length_scale[2]))

}

compute_kernel_tilde<-function(a)
{
x<-a[1:2]
y<-a[3:4]
w1<- a[5]
w2<-a[6]


d<- sum((x-y)^2)
u11<- (w1_max-w1)*sqrt(2/length_scale[1])
l11<- (w1_min-w1)*sqrt(2/length_scale[1])
u21<- (w2_max-w2)*sqrt(2/length_scale[1])
l21<- (w2_min-w2)*sqrt(2/length_scale[1])

u12<- (w1_max-w1)*sqrt(2/length_scale[2])
l12<- (w1_min-w1)*sqrt(2/length_scale[2])
u22<- (w2_max-w2)*sqrt(2/length_scale[2])
l22<- (w2_min-w2)*sqrt(2/length_scale[2])

c(sigmaf[1]*exp(-d/length_scale[1])*(  sqrt(length_scale[1]*pi)/(2*w1_max) *( pnorm(u11)-pnorm(l11))   ) *(  sqrt(length_scale[1]*pi)/(2*w2_max) *( pnorm(u21)-pnorm(l21))   ) ,
sigmaf[2]*exp(-d/length_scale[2])*(  sqrt(length_scale[2]*pi)/(2*w1_max) *( pnorm(u12)-pnorm(l12))   ) *(  sqrt(length_scale[2]*pi)/(2*w2_max) *( pnorm(u22)-pnorm(l22))   ) )
}

compute_kernel_tilde2<-function(a)
{
x<-a[1:2]
y<-a[3:4]




d<- sum((x-y)^2)
A11<- (1/(4*w1_max^2))*sqrt(length_scale[1]/2)*sqrt(length_scale[1]*pi)
A21<- (1/(4*w2_max^2))*sqrt(length_scale[1]/2)*sqrt(length_scale[1]*pi)
A12<- (1/(4*w1_max^2))*sqrt(length_scale[2]/2)*sqrt(length_scale[2]*pi)
A22<- (1/(4*w2_max^2))*sqrt(length_scale[2]/2)*sqrt(length_scale[2]*pi)

B11<- 2*w1_max*sqrt(2/length_scale[1])
B21<- 2*w2_max*sqrt(2/length_scale[1])
B12<- 2*w1_max*sqrt(2/length_scale[2])
B22<- 2*w2_max*sqrt(2/length_scale[2])

c(sigmaf[1]*exp(-d/length_scale[1])* A11*(B11*(pnorm(B11)-pnorm(-B11))+2*standard_normal_pdf(B11)-2*standard_normal_pdf(0))*A21*(B21*(pnorm(B21)-pnorm(-B21))+2*standard_normal_pdf(B21)-2*standard_normal_pdf(0)),
sigmaf[2]*exp(-d/length_scale[2])* A12*(B12*(pnorm(B12)-pnorm(-B12))+2*standard_normal_pdf(B12)-2*standard_normal_pdf(0))*A22*(B22*(pnorm(B22)-pnorm(-B22))+2*standard_normal_pdf(B22)-2*standard_normal_pdf(0)))


}



compute_ParEGO_kernel<-function(a)
{
x1<-a[1]
y1<-a[2]
x2<-a[3]
y2<-a[4]

ParEGO_sigmaf*exp(-((x1-x2)^2+(y1-y2)^2)/ParEGO_length)
}

generate_ParEGO_weight_vector<-function(od)
{
A<-matrix(runif(10000*(od-1)),10000,(od-1))
	while(sum(rowSums(A)<=1) < 0.5)
	{
	A<-matrix(runif(10000*(od-1)),10000,(od-1))
	}

a<-A[c(1:10000)[rowSums(A)<=1][1],]
c(a,1-sum(a))
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
set.seed(1)
results <- nsga2R(fn=negative_true_Bayes_risk_function, varNo=2, objDim=2, lowerBounds=c(x1_min,x2_min), upperBounds=c(x1_max,x2_max),
popSize=500, tourSize=2, generations=100, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)
#plot(results$objectives)






true_PF<- -(results$objectives)

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
#X_num<-c()

simple_PHV_regret<-c()
simple_e_accuracy<-c()
simple_e_accuracy2<-c()


inference_PHV_regret<-c()
inference_e_accuracy<-c()
inference_e_accuracy2<-c()


I_1<- matrix(numeric(iteration_num^2),iteration_num,iteration_num)
I_2<- matrix(numeric(iteration_num^2),iteration_num,iteration_num)


new_input<- runif(2,c(x1_min,x2_min),c(x1_max,x2_max))
new_w1<- runif(max(naive_sample_num,1000),w1_min,w1_max)
new_w2<- runif(max(naive_sample_num,1000),w2_min,w2_max)

new_Y_all<- t(apply(cbind(new_input[1],new_input[2],new_w1,new_w2),1,true_function))+mvrnorm(max(naive_sample_num,1000),c(0,0),diag(sigman))

new_Y<- colMeans(new_Y_all[1:naive_sample_num,])
estimated_var_1<- var(new_Y_all[,1])
estimated_var_2<- var(new_Y_all[,2])


K1<- matrix(numeric(iteration_num^2),iteration_num,iteration_num)	
K2<- matrix(numeric(iteration_num^2),iteration_num,iteration_num)	

	

	for(ww in 1:iteration_num)
	{
	I_1[ww,ww]<- (estimated_var_1+sigman[1])/naive_sample_num
	I_2[ww,ww]<- (estimated_var_2+sigman[2])/naive_sample_num

	input<-rbind(input,new_input)	
	#X_num<-c(X_num,new_X_num)
	Y<-rbind(Y,new_Y)

	k_vec<-apply(cbind(new_input[1],new_input[2],input),1,compute_kernel_naive)	

	K1[ww,1:ww]<- k_vec[1,]
	K2[ww,1:ww]<- k_vec[2,]

		if(ww==1)
		{
		G1<- K1[1:ww,1:ww]
		G2<- K2[1:ww,1:ww]
		} else {
			G1<- K1[1:ww,1:ww]+t(K1[1:ww,1:ww])-diag(diag(K1[1:ww,1:ww]))
			G2<- K2[1:ww,1:ww]+t(K2[1:ww,1:ww])-diag(diag(K2[1:ww,1:ww]))
		}
	C1<- G1+ I_1[1:ww,1:ww]
	C2<- G2+ I_2[1:ww,1:ww]


	C1_inv<- solve(C1)
	C2_inv<- solve(C2)

	
	C1_inv_Y1<- as.vector(C1_inv%*%Y[,1])
	C2_inv_Y2<- as.vector(C2_inv%*%Y[,2])

	post_mean_function<-function(x)
	{
	A<- cbind(x[1],x[2],input)
	B<- rbind(t(apply(A,1,compute_kernel_naive)))
	C<- cbind(C1_inv_Y1,C2_inv_Y2)
	colSums(B*C)
	}

	post_var_function<-function(x)
	{
	A<- cbind(x[1],x[2],input)
	B<- apply(A,1,compute_kernel_naive)
	c(as.vector(abs(sigmaf[1]-t(B[1,])%*%C1_inv%*%B[1,])	),as.vector(abs(sigmaf[2]-t(B[2,])%*%C2_inv%*%B[2,])	))
	}

	ucb_function<-function(x)
	{
	post_mean_function(x)+beta*sqrt(post_var_function(x))
	}

	lcb_function<-function(x)
	{
	post_mean_function(x)-beta*sqrt(post_var_function(x))
	}

	negative_lcb_function<-function(x)
	{
	-lcb_function(x)
	}

	
	predicted_PFn<-   compute_PF(t(apply(rbind(input[,1:2]),1,true_Bayes_risk_function)))
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

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)

		new_input<- runif(2,c(x1_min,x2_min),c(x1_max,x2_max))

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		
		}


		if(AF_num==2)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		aaa<- cbind(runif(100,x1_min,x1_max),runif(100,x2_min,x2_max))
		bbb<- colSums(apply(aaa,1,post_var_function))
		ini<- aaa[order(bbb,decreasing=T)[1],]
				
		AF_x<- function(a)
		{
		-sum(post_var_function(a))
		}

		resu<- optim(ini,AF_x,method = c("L-BFGS-B"),lower = c(x1_min,x2_min), upper = c(x1_max,x2_max))
		new_input<- resu$par

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		
		}


		if(AF_num==3)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		current_PF<- rbind(rbind(compute_PF(Y))[,1:obj_dim])
		current_PF_cell<-create_cells(current_PF,apply(rbind(current_PF,ref_point),2,min),apply(current_PF,2,max),c(1:length(current_PF[,1])))
		current_PHV<- sum(apply(rbind(current_PF_cell[,(obj_dim+1):(2*obj_dim)]-current_PF_cell[,1:(obj_dim)]),1,prod))
		EHVI_sample<- mvrnorm(EHVI_sample_num,numeric(obj_dim),diag(obj_dim))

		AF_x<- function(a)
		{
		mu<- post_mean_function(a)
		sd<- sqrt(post_var_function(a))
		
		y1<- mu[1]+sd[1]*EHVI_sample[,1]
		y2<- mu[2]+sd[2]*EHVI_sample[,2]
		newY<- cbind(y1,y2)
		resu<-numeric(EHVI_sample_num)
			for(i in 1:EHVI_sample_num)
			{
			add_Y_cPF<- rbind(current_PF,newY[i,])
			new_cPF<- rbind(rbind(compute_PF(add_Y_cPF))[,1:obj_dim])
			new_cPF_cell<- create_cells(new_cPF,apply(rbind(new_cPF,ref_point),2,min),apply(new_cPF,2,max),c(1:length(new_cPF[,1])))
			new_cPF_PHV<- sum(apply(rbind(new_cPF_cell[,(obj_dim+1):(2*obj_dim)]-new_cPF_cell[,1:(obj_dim)]),1,prod))
			resu[i]<- new_cPF_PHV-current_PHV
			}
		-mean(resu)
		}


		aaa<- cbind(runif(100,x1_min,x1_max),runif(100,x2_min,x2_max))
		bbb<- apply(aaa,1,AF_x)
		ini<- aaa[order(bbb)[1],]
				
		

		resu<- optim(ini,AF_x,method = c("L-BFGS-B"),lower = c(x1_min,x2_min), upper = c(x1_max,x2_max))
		new_input<- resu$par

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		
		}


		if(AF_num==4)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		results <- nsga2R(fn=negative_lcb_function, varNo=2, objDim=2, lowerBounds=c(x1_min,x2_min), upperBounds=c(x1_max,x2_max),
		popSize=50, tourSize=2, generations=50, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)

		lcb_PF<- -(results$objectives)



		AF_x<- function(a)
		{
		b<- ucb_function	(a)
		bb<- lcb_function	(a)

		c<-(min(apply(b-t(lcb_PF),2,max))>0)*1
		(-1)*sum((b-bb)^2)*c
		
		}


		aaa<- cbind(runif(100,x1_min,x1_max),runif(100,x2_min,x2_max))
		bbb<- apply(aaa,1,AF_x)
		ini<- aaa[order(bbb)[1],]
				
		

		resu<- optim(ini,AF_x,method = c("L-BFGS-B"),lower = c(x1_min,x2_min), upper = c(x1_max,x2_max))
		new_input<- resu$par

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		
		}


		if(AF_num==5)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		results <- nsga2R(fn=negative_lcb_function, varNo=2, objDim=2, lowerBounds=c(x1_min,x2_min), upperBounds=c(x1_max,x2_max),
		popSize=50, tourSize=2, generations=50, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)

		lcb_PF<- -(results$objectives)



		AF_x<- function(a)
		{
		b<- ucb_function	(a)


		
		-min(apply(b-t(lcb_PF),2,max))
		
		}


		aaa<- cbind(runif(100,x1_min,x1_max),runif(100,x2_min,x2_max))
		bbb<- apply(aaa,1,AF_x)
		ini<- aaa[order(bbb)[1],]
				
		

		resu<- optim(ini,AF_x,method = c("L-BFGS-B"),lower = c(x1_min,x2_min), upper = c(x1_max,x2_max))
		new_input<- resu$par

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		
		}



		if(AF_num==6)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		current_PF<- rbind(rbind(compute_PF(Y))[,1:obj_dim])
		current_PF_cell<-create_cells(current_PF,apply(rbind(current_PF,ref_point),2,min),apply(current_PF,2,max),c(1:length(current_PF[,1])))
		current_PHV<- sum(apply(rbind(current_PF_cell[,(obj_dim+1):(2*obj_dim)]-current_PF_cell[,1:(obj_dim)]),1,prod))
		EMMDI_sample<- mvrnorm(EMMDI_sample_num,numeric(obj_dim),diag(obj_dim))

		AF_x<- function(a)
		{
		mu<- post_mean_function(a)
		sd<- sqrt(post_var_function(a))
		
		y1<- mu[1]+sd[1]*EMMDI_sample[,1]
		y2<- mu[2]+sd[2]*EMMDI_sample[,2]
		newY<- cbind(y1,y2)
		 -compute_EMMDI(newY,current_PF)
		}


		aaa<- cbind(runif(100,x1_min,x1_max),runif(100,x2_min,x2_max))
		bbb<- apply(aaa,1,AF_x)
		ini<- aaa[order(bbb)[1],]
				
		

		resu<- optim(ini,AF_x,method = c("L-BFGS-B"),lower = c(x1_min,x2_min), upper = c(x1_max,x2_max))
		new_input<- resu$par

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

		
		}



	######################################################################
	results <- nsga2R(fn=negative_lcb_function, varNo=2, objDim=2, lowerBounds=c(x1_min,x2_min), upperBounds=c(x1_max,x2_max),
		popSize=50, tourSize=2, generations=50, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)

		reported_X<-results$parameters
		reported_Y<-rbind(t(apply(reported_X,1,true_Bayes_risk_function)))

		inference_predicted_PF_PFn<-compute_PF(rbind(t(apply(reported_X,1,true_Bayes_risk_function))))
		inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])
		
		inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
		inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


		inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
		inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

		inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
		inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

		inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,reported_Y),0))
		inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]


	set.seed(((seed-1)*total_seed_num+w)*1000+ww)

	new_w1<- runif(max(naive_sample_num,1000),w1_min,w1_max)
	new_w2<- runif(max(naive_sample_num,1000),w2_min,w2_max)

	new_Y_all<- t(apply(cbind(new_input[1],new_input[2],new_w1,new_w2),1,true_function))+mvrnorm(max(naive_sample_num,1000),c(0,0),diag(sigman))

	new_Y<- colMeans(new_Y_all[1:naive_sample_num,])
	estimated_var_1<- var(new_Y_all[,1])
	estimated_var_2<- var(new_Y_all[,2])


		
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




