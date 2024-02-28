library(MASS)
library(lpSolve)

seed<-1
total_seed_num<-100
iteration_num<-300
PHV_sample_num<-10000
AF_num<-14			#1:Random, 2:US, 3:EHVI, 4:e-PAL, 5:Proposed, 6:modified_proposed, 7:ParEGO, 8:SMSego, 9:MESMO, 10:PFES, 11:AD-GP-UCB, 12:MOBO_RS_UCB, 13:MOBO_RS_TS, 14:EMMDI


obj_num<-c(1,2)
obj_dim<-length(obj_num)
input_dim<-2
grid_num<-50
total_input_num<-grid_num^input_dim

lower<- -5
upper<- 5

sigmaf<-c(2,2,2,2,2)
sigman<-10^(-6)
length_scale<-c(2,2,2,2,2)
beta<-3


ParEGO_rho<-0.05
ParEGO_length<-2
ParEGO_sigmaf<-2
ParEGO_sigman<-10^(-8)

SMSego_sample_num<-50

EHVI_sample_num<-20

MESMO_sample_num<-10
RFM_sample_num<-500

PFES_sample_num<-10
max_PF_num<-50


EMMDI_sample_num<-20


###############################################################
##################     True function    #######################
###############################################################

true_function<-function(x,y,l)
{
	if(l==1)#Minus Booth function
	{
	f= (-1*( (x+2*y-7)^2+(2*x+y-5)^2 )-(-157.35))/sqrt(28896.11)
	}

	if(l==2)#Minus Matyas function
	{
	f= (-1*( 0.26*(x^2+y^2)-0.48*x*y )-(-4.3342))/sqrt(23.52052)
	}

	if(l==3)#Minus Himmelblau's function
	{
	f= (-1*( (x^2+y-11)^2+(x+y^2-7)^2 )-(-136.71))/sqrt(12503.63)
	}

	if(l==4)#Minus McCormic function
	{
	f= (-1*( sin(x+y)+(x-y)^2-1.5*x+2.5*y+1 )-(-17.67))/sqrt(460.573)
	}

	if(l==5)#Minus Styblinski-Tang function
	{
	f= (-1*( (x^4-16*x^2+5*x)/2+(y^4-16*y^2+5*y)/2 )-(8.309994))/sqrt(2059.096)
	}
f
}

###############################################################
###############################################################
###############################################################



###############################################################
##################   Input and Output   #######################
###############################################################

x1_all<- seq(lower,upper,length=grid_num)
x2_all<- seq(lower,upper,length=grid_num)
all_input<-cbind((numeric(grid_num)+1)%x%x1_all,x2_all%x%(numeric(grid_num)+1))

true_function_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
for(i in 1:obj_dim)
{
true_function_mat[,i]<- true_function(all_input[,1],all_input[,2],obj_num[i])
}

ref_point<-apply(true_function_mat,2,min)
max_point<-apply(true_function_mat,2,max)
vol_refp_maxp<- prod(max_point-ref_point)

###############################################################
###############################################################
###############################################################



###############################################################
##################       Function       #######################
###############################################################

compute_kernel<-function(a)
{
x1<-a[1]
y1<-a[2]
x2<-a[3]
y2<-a[4]
l<-a[5]

sigmaf[l]*exp(-((x1-x2)^2+(y1-y2)^2)/length_scale[l])
}

compute_kernel_deri1<-function(a)
{
x1<-a[1]
y1<-a[2]
x2<-a[3]
y2<-a[4]
l<-a[5]

sigmaf[l]*exp(-((x1-x2)^2+(y1-y2)^2)/length_scale[l])*(  -2*(x1-x2)/length_scale[l]    )
}


compute_kernel_deri2<-function(a)
{
x1<-a[1]
y1<-a[2]
x2<-a[3]
y2<-a[4]
l<-a[5]

sigmaf[l]*exp(-((x1-x2)^2+(y1-y2)^2)/length_scale[l])*(  -2*(y1-y2)/length_scale[l]    )
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

standard_normal_pdf<-function(x)
{
(2*pi)^(-1/2)*exp(-x^2/2)
}


#compute_PF<-function(A)
#{
#	if(length(A[,1])==1)
#	{
#	B<-cbind(A,1)
#	} else {
#	l<- length(A[,1])
#	B<-c()
#		for(j in 1:l)
#		{
#		a<-A[j,]
#		AA<-A[-j,]
#		AAA<-t(AA)-(a)
#		aa<-max(apply(AAA,2,min))
#			if(aa<0)
#			{
#			B<-rbind(B,c(a,j))
#			}
#		}
#
#	}
#B
#}

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


compute_PF_cell_dim2<-function(PF)
{
PF<-rbind(PF)
rp<-apply(cbind(apply(PF,2,min),c(-10^(10),-10^(10))),1,min)
num<- order(PF[,2],decreasing=T)
PF<-rbind(PF[num,])
cbind(c(rp[1],PF[,1][-length(num)]),rp[2],PF)
}

compute_PF_cell_dim3<-function(PF)
{
PF<-rbind(PF)
rp<-apply(cbind(apply(PF,2,min),c(-10^(10),-10^(10),-10^(10))),1,min)
num<- order(PF[,3],decreasing=T)
PF<-rbind(PF[num,])
A<-c()
	for(i in 1:length(num))
	{
	PF2<- rbind(PF[,1:2])[1:i,]
	PF2<-rbind(PF2)
	PF2<-rbind(compute_PF(PF2)[,1:2])
	B<-compute_PF_cell_dim2(PF2)
		if(i==length(num))
		{
		B<- cbind(rbind(B[,1:2]),rp[3],rbind(B[,3:4]),PF[,3][i])
		} else {
		B<- cbind(rbind(B[,1:2]),PF[,3][(i+1)],rbind(B[,3:4]),PF[,3][i])
		}
	A<-rbind(A,B)
	}
A
}

compute_PF_cell_dim4<-function(PF)
{
PF<-rbind(PF)
rp<-apply(cbind(apply(PF,2,min),c(-10^(10),-10^(10),-10^(10),-10^(10))),1,min)
num<- order(PF[,4],decreasing=T)
PF<-rbind(PF[num,])
A<-c()
	for(i in 1:length(num))
	{
	PF3<- rbind(PF[,1:3])[1:i,]
	PF3<-rbind(PF3)
	PF3<-rbind(compute_PF(PF3)[,1:3])
	B<-compute_PF_cell_dim3(PF3)
		if(i==length(num))
		{
		B<- cbind(rbind(B[,1:3]),rp[4],rbind(B[,4:6]),PF[,4][i])
		} else {
		B<- cbind(rbind(B[,1:3]),PF[,4][(i+1)],rbind(B[,4:6]),PF[,4][i])
		}
	A<-rbind(A,B)
	}
A
}

compute_PF_cell_dim5<-function(PF)
{
PF<-rbind(PF)
rp<-apply(cbind(apply(PF,2,min),c(-10^(10),-10^(10),-10^(10),-10^(10),-10^(10))),1,min)
num<- order(PF[,5],decreasing=T)
PF<-rbind(PF[num,])
A<-c()
	for(i in 1:length(num))
	{
	PF4<- rbind(PF[,1:4])[1:i,]
	PF4<-rbind(PF4)
	PF4<-rbind(compute_PF(PF4)[,1:4])
	B<-compute_PF_cell_dim4(PF4)
		if(i==length(num))
		{
		B<- cbind(rbind(B[,1:4]),rp[5],rbind(B[,5:8]),PF[,5][i])
		} else {
		B<- cbind(rbind(B[,1:4]),PF[,5][(i+1)],rbind(B[,5:8]),PF[,5][i])
		}
	A<-rbind(A,B)
	}
A
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
#A<-matrix(numeric(2*l*M),M,(2*l))
#alpha<-matrix(numeric(2*l*M),M,(2*l))
#UP<-matrix(numeric(l*M),M,(l))
#LO<-matrix(numeric(l*M),M,(l))
#Zml<-matrix(numeric(l*M),M,(l))
#Zm<-numeric(M)
#Gammaml<-matrix(numeric(l*M),M,(l))


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


#compute_PFES2<-function(pm,pv,cell)
#{
#M<-length(cell[,1])
#l<-length(cell[1,])/2
#Z<-0
#ZmGamma<-0
#A<-matrix(numeric(2*l*M),M,(2*l))
#alpha<-matrix(numeric(2*l*M),M,(2*l))
#UP<-matrix(numeric(l*M),M,(l))
#LO<-matrix(numeric(l*M),M,(l))
#Zml<-matrix(numeric(l*M),M,(l))
#Zm<-numeric(M)
#Gammaml<-matrix(numeric(l*M),M,(l))
#
#
#A<-t(c(pm,pm))%x%(numeric(M)+1)
#alpha<- (cell-A)/sqrt(pv)
#UP<- rbind(alpha[,(l+1):(2*l)])
#LO<-rbind(alpha[,1:l])
#
#Zml<-(pnorm(UP)-pnorm(LO))
#Zml2<-Zml
#for(i in 1:l)
#{
#Zml2[,l]<- apply(cbind(Zml[,-l]),1,prod)
#}
#
#Zm<- apply(Zml,1,prod)
#Z<-sum(Zm)
#Gammaml<- (LO*standard_normal_pdf(LO)-UP*standard_normal_pdf(UP))
#
#ZmGamma<- sum(rowSums(Zml2*Gammaml))/2
#
#if(Z>1)
#{
#Z<-1
#}
#
#   -log(Z) - ZmGamma/Z
#
#}


score_function<-function(x)
{
(1+exp(-x))^(-1)
}

deri_score_function<-function(x)
{
exp(-x)*(1+exp(-x))^(-2)
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




###############################################################
###############################################################
###############################################################


###############################################################
##################     Pareto Front     #######################
###############################################################

true_PF_PFn<-compute_PF(true_function_mat)
true_PF<-true_PF_PFn[,1:obj_dim]
true_PFn<-true_PF_PFn[,(obj_dim+1)]

#set.seed(1)
true_PF_cell<- create_cells(true_PF,ref_point,max_point,c(1:length(true_PF[,1])))
true_PHV<- sum(apply(rbind(rbind(true_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(true_PF_cell)[,1:(obj_dim)]),1,prod))

###############################################################
###############################################################
###############################################################


###############################################################
##################     Random feature       ###################
###############################################################
set.seed(1)
RFM_mat<-  matrix(numeric(total_input_num*RFM_sample_num),total_input_num,RFM_sample_num)%x%numeric(obj_dim)
for(i in 1:obj_dim)
{
RFM_w_mat<-  mvrnorm(RFM_sample_num,numeric(input_dim),(1/(2/length_scale[obj_num[i]]))*diag(input_dim))
RFM_b_vec<- runif(RFM_sample_num,0,2*pi)
	for(j in 1:total_input_num)
	{
	RFM_mat[((i-1)*total_input_num+j),]<-     sqrt(2*sigmaf[obj_num[i]]/RFM_sample_num)*cos(as.vector(RFM_w_mat%*%all_input[j,])+RFM_b_vec)
	}
}



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

X_num_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
AF_computation_time_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)

###############################################################
###############################################################
###############################################################


###############################################################
##################     Optimization     #######################
###############################################################

Kvec_mat<-matrix(numeric(total_input_num*iteration_num),total_input_num,iteration_num)%x%numeric(obj_dim)
post_mean_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
post_var_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
ucb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
lcb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
post_cov_mat<-diag(total_input_num)%x%numeric(obj_dim)

svd_mat<-diag(total_input_num)%x%numeric(obj_dim)




for(w in 1:total_seed_num)
{
set.seed((seed-1)*total_seed_num+w)

X<-c()
Y<-c()
X_num<-c()

simple_PHV_regret<-c()
simple_e_accuracy<-c()
simple_e_accuracy2<-c()


inference_PHV_regret<-c()
inference_e_accuracy<-c()
inference_e_accuracy2<-c()


new_X_num<- sample(c(1:total_input_num),1)
new_X<- all_input[new_X_num,]
new_Y<- true_function_mat[new_X_num,]+rnorm(obj_dim,0,sqrt(sigman))
	
	for(ww in 1:iteration_num)
	{
	X<-rbind(X,new_X)	
	X_num<-c(X_num,new_X_num)
	Y<-rbind(Y,new_Y)

		for(www in 1:obj_dim)
		{
		Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),ww]<-t(matrix(apply(cbind(all_input,t(new_X)%x%(numeric(total_input_num)+1),www),1,compute_kernel),1,total_input_num))

		partial_Kvec_mat<-cbind(Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),1:ww])
		Cinv<-solve(partial_Kvec_mat[X_num,]+sigman*diag(ww))

		post_mean_mat[,www]<-as.vector(partial_Kvec_mat%*%Cinv%*%Y[,www])
		post_var_mat[,www]<- abs(sigmaf[www]-rowSums((partial_Kvec_mat%*%Cinv)*partial_Kvec_mat))

		}

	ucb_mat<- post_mean_mat+beta*sqrt(post_var_mat)
	lcb_mat<- post_mean_mat-beta*sqrt(post_var_mat)

	predicted_PF_PFn<-compute_PF(rbind(true_function_mat[X_num,]))
	predicted_PF<-rbind(predicted_PF_PFn[,1:obj_dim])
	predicted_PFn<-X_num[predicted_PF_PFn[,(obj_dim+1)]]


	if(obj_dim==5)
	{
	simple_PHV_regret<- c(simple_PHV_regret,max(true_PHV-compute_PHV(predicted_PF),0))
	simple_PHV_regret_mat[w,ww]<-	min(simple_PHV_regret)

	simple_e_accuracy<- c(simple_e_accuracy,compute_e_accuracy(true_PF,predicted_PF))
	simple_e_accuracy_mat[w,ww]<-	min(simple_e_accuracy)

	simple_e_accuracy2<- c(simple_e_accuracy2,compute_e_accuracy2(true_PF,predicted_PF))
	simple_e_accuracy2_mat[w,ww]<-	min(simple_e_accuracy2)


	lcb_PF_num<-compute_PF(lcb_mat)[,(obj_dim+1)]
	reported_Y<-rbind(true_function_mat[lcb_PF_num,])
	inference_predicted_PF_PFn<-compute_PF(rbind(true_function_mat[lcb_PF_num,]))
	inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])
	inference_predicted_PFn<-lcb_PF_num[inference_predicted_PF_PFn[,(obj_dim+1)]]


	inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-compute_PHV(inference_predicted_PF),0))
	inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

	inference_e_accuracy<- c(inference_e_accuracy,compute_e_accuracy(true_PF,inference_predicted_PF))
	inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

	inference_e_accuracy2<- c(inference_e_accuracy2,compute_e_accuracy2(true_PF,reported_Y))
	inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]
	} else {
	predicted_PF_cell<-create_cells(predicted_PF,ref_point,apply(predicted_PF,2,max),c(1:length(predicted_PF[,1])))

	predicted_PF_PHV<-  sum(apply(rbind(predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-predicted_PF_cell[,(1):(obj_dim)]),1,prod))
	
	
	simple_PHV_regret<- c(simple_PHV_regret,max(true_PHV-predicted_PF_PHV,0))
	simple_PHV_regret_mat[w,ww]<-	min(simple_PHV_regret)

	simple_e_accuracy<- c(simple_e_accuracy,compute_e_accuracy(true_PF,predicted_PF))
	simple_e_accuracy_mat[w,ww]<-	min(simple_e_accuracy)

	simple_e_accuracy2<- c(simple_e_accuracy2,compute_e_accuracy2(true_PF,predicted_PF))
	simple_e_accuracy2_mat[w,ww]<-	min(simple_e_accuracy2)


	lcb_PF_num<-compute_PF(lcb_mat)[,(obj_dim+1)]
	reported_Y<-rbind(true_function_mat[lcb_PF_num,])
	inference_predicted_PF_PFn<-compute_PF(rbind(true_function_mat[lcb_PF_num,]))
	inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])
	inference_predicted_PFn<-lcb_PF_num[inference_predicted_PF_PFn[,(obj_dim+1)]]

	inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
	inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


	inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
	inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

	inference_e_accuracy<- c(inference_e_accuracy,compute_e_accuracy(true_PF,inference_predicted_PF))
	inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

	inference_e_accuracy2<- c(inference_e_accuracy2,compute_e_accuracy2(true_PF,reported_Y))
	inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]


	}

	X_num_mat[w,ww]<-new_X_num

	################## Acquisition function ##############################
		if(AF_num==1)
		{
		starttime<-Sys.time()

		new_X_num<- sample(c(1:total_input_num),1)

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==2)
		{
		starttime<-Sys.time()

		AF<- rowSums(post_var_mat)
		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}

		if(AF_num==3)
		{
		starttime<-Sys.time()
		estimated_PF<- rbind((compute_PF(Y))[,1:obj_dim])
		AF<-numeric(total_input_num)
			for(abc in 1:total_input_num)
			{
				for(abcabc in 1:EHVI_sample_num)
				{
				sample_f<-mvrnorm(1,post_mean_mat[abc,],diag(post_var_mat[abc,]))
				AF[abc]<-AF[abc]+	compute_SMSego_PHVI(sample_f,estimated_PF,EHVI_sample_num)
				}
			}		
		AF<-AF/EHVI_sample_num
		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}




		if(AF_num==4)
		{
		starttime<-Sys.time()

		lcb_PF_num<-compute_PF(lcb_mat)[,(obj_dim+1)]
		AF1<-compute_dist(ucb_mat,rbind(lcb_mat[lcb_PF_num,]))
		AF2<-rowSums((ucb_mat-lcb_mat)^2)*(AF1>=0)

		new_X_num<- order(AF2,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==5)
		{
		starttime<-Sys.time()

		lcb_PF_num<-compute_PF(lcb_mat)[,(obj_dim+1)]
		AF1<-compute_dist(ucb_mat,rbind(lcb_mat[lcb_PF_num,]))

		new_X_num<- order(AF1,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==6)
		{
		starttime<-Sys.time()

		lcb_PF_num<-compute_PF(rbind(lcb_mat[X_num,]))[,(obj_dim+1)]
		AF1<-compute_dist(ucb_mat,rbind(lcb_mat[X_num[lcb_PF_num],]))

		new_X_num<- order(AF1,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==7)
		{
		starttime<-Sys.time()

		ParEGO_weight_vector<- generate_ParEGO_weight_vector(obj_dim)
		ParEGO_Y<-  apply(Y*(numeric(ww)+1)%x%t(ParEGO_weight_vector),1,max)+as.vector(ParEGO_rho*Y%*%ParEGO_weight_vector)
		ParEGO_current_best<-max(ParEGO_Y)
		
		ParEGO_kvec<- t(matrix(apply(cbind(all_input%x%(numeric(ww)+1),(numeric(total_input_num)+1)%x%X),1,compute_ParEGO_kernel),ww,total_input_num))
		ParEGO_Kmat<- ParEGO_kvec[X_num,]	
		ParEGO_Cinv<- solve(ParEGO_Kmat+ParEGO_sigman*diag(ww))	

		ParEGO_post_mean<- as.vector(ParEGO_kvec%*%ParEGO_Cinv%*%ParEGO_Y)
		ParEGO_post_var<- abs(ParEGO_sigmaf-rowSums((ParEGO_kvec%*%ParEGO_Cinv)*ParEGO_kvec))

		ParEGO_Z<- (ParEGO_post_mean-ParEGO_current_best)/sqrt(ParEGO_post_var)
		AF<-sqrt(ParEGO_post_var)*(ParEGO_Z*pnorm(ParEGO_Z)+standard_normal_pdf(ParEGO_Z))

		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==8)
		{
		starttime<-Sys.time()
		estimated_PF<- rbind((compute_PF(Y))[,1:obj_dim])
		AF<-numeric(total_input_num)
			for(abc in 1:total_input_num)
			{
			AF[abc]<-	compute_SMSego_PHVI(ucb_mat[abc,],estimated_PF,SMSego_sample_num)
			}		

		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}

		
		if(AF_num==9)
		{
		starttime<-Sys.time()
		sample_y_mat<-matrix(numeric(total_input_num*MESMO_sample_num*obj_dim),total_input_num*obj_dim,MESMO_sample_num)
		#RFM_post_mean_mat<-post_mean_mat*0 
		#RFM_post_var_mat<- post_var_mat*0 

			for(www in 1:obj_dim)
			{
			RFM_X<- RFM_mat[((www-1)*total_input_num+1):(www*total_input_num),]
			current_RFM_X<-rbind(RFM_X[X_num,])
			inv_tXX_I<- solve(t(current_RFM_X)%*%current_RFM_X+sigman*diag(RFM_sample_num))
			post_coefficient<-inv_tXX_I%*%t(current_RFM_X)%*%Y[,www]
			post_RFM_cov<- sigman*inv_tXX_I
			#svd_u_d<-svd(post_RFM_cov)
			#generated_coeficient<-as.vector(post_coefficient)+(svd_u_d$u)%*% diag(sqrt(abs(svd_u_d$d)))%*%matrix(rnorm(RFM_sample_num*MESMO_sample_num,0,1),RFM_sample_num,MESMO_sample_num)
			generated_coeficient<-t(mvrnorm(MESMO_sample_num,post_coefficient,(post_RFM_cov+10^(-8)*diag(RFM_sample_num))))
			sample_y_mat[((www-1)*total_input_num+1):(www*total_input_num),]<-RFM_X%*%generated_coeficient
			#RFM_post_mean_mat[,www]<-		as.vector(RFM_X%*%post_coefficient)
			#RFM_post_var_mat[,www]<-	 abs(rowSums((RFM_X%*%post_RFM_cov)*RFM_X) )

			}
		
			AF<-numeric(total_input_num)

			for(sa in 1:MESMO_sample_num)
			{
			sample_PF<-rbind(compute_PF(matrix(sample_y_mat[,sa],total_input_num,obj_dim))[,1:obj_dim])
			
				for(id in 1:total_input_num)
				{
				AF[id]<- AF[id]+compute_MESMO(sample_PF,post_mean_mat[id,],post_var_mat[id,])
				}

			}

		if(max(AF)==0)
		{
		new_X_num<- sample(c(1:total_input_num),1)
		} else {
		new_X_num<- order(AF,decreasing=T)[1]
		}

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}



		if(AF_num==10)
		{
		starttime<-Sys.time()
		sample_y_mat<-matrix(numeric(total_input_num*PFES_sample_num*obj_dim),total_input_num*obj_dim,PFES_sample_num)
			#RFM_post_mean_mat<-post_mean_mat*0 
			#RFM_post_var_mat<- post_var_mat*0 

			for(www in 1:obj_dim)
			{
			RFM_X<- RFM_mat[((www-1)*total_input_num+1):(www*total_input_num),]
			current_RFM_X<-rbind(RFM_X[X_num,])
			inv_tXX_I<- solve(t(current_RFM_X)%*%current_RFM_X+sigman*diag(RFM_sample_num))
			post_coefficient<-inv_tXX_I%*%t(current_RFM_X)%*%Y[,www]
			post_RFM_cov<- sigman*inv_tXX_I
			#svd_u_d<-svd(post_RFM_cov)
			#generated_coeficient<-as.vector(post_coefficient)+(svd_u_d$u)%*% diag(sqrt(abs(svd_u_d$d)))%*%matrix(rnorm(RFM_sample_num*MESMO_sample_num,0,1),RFM_sample_num,MESMO_sample_num)
			generated_coeficient<-t(rbind(mvrnorm(PFES_sample_num,post_coefficient,(post_RFM_cov+10^(-8)*diag(RFM_sample_num)))))
			sample_y_mat[((www-1)*total_input_num+1):(www*total_input_num),]<-RFM_X%*%generated_coeficient
				#RFM_post_mean_mat[,www]<-		as.vector(RFM_X%*%post_coefficient)
				#RFM_post_var_mat[,www]<-	 abs(rowSums((RFM_X%*%post_RFM_cov)*RFM_X) )
			}
		
			AF<-numeric(total_input_num)
			#AF2<-numeric(total_input_num)


			for(sa in 1:PFES_sample_num)
			{
			SAMPLE_Y_MAT<- matrix(sample_y_mat[,sa],total_input_num,obj_dim)

			sample_PF<-rbind(compute_PF(matrix(sample_y_mat[,sa],total_input_num,obj_dim))[,1:obj_dim])
			dist_id<- c(1:length(sample_PF[,1]))	
			
				if(length(sample_PF[,1])>max_PF_num && obj_dim>3)
				{
				add_sample_PF<- rbind(sample_PF,apply(sample_PF,2,min))
				dist_mat<- matrix(rowSums(((add_sample_PF[,1]*0+1)%x%add_sample_PF-add_sample_PF%x%(add_sample_PF[,1]*0+1))^2),length(add_sample_PF[,1]),length(add_sample_PF[,1]))
				dist_id<-c(length(add_sample_PF[,1]))	
					for(wwww in 1:min(max_PF_num,length(sample_PF[,1])))
					{
					dist_id<- c(dist_id,order(apply(cbind(dist_mat[,dist_id]),1,min),decreasing=T)[1])
					}
				dist_id<-dist_id[-1]

					for(wwwww in 1:length(X_num))
					{
					S<-apply(t(sample_PF[dist_id,])-SAMPLE_Y_MAT[X_num[wwwww],],2,min)
						if(max(S)<0)
						{
						dist_id<-c(dist_id,order(apply(t(sample_PF)-SAMPLE_Y_MAT[X_num[wwwww],],2,min),decreasing=T)[1])
						}
					}

				}



			#cell<- create_cells(sample_PF,(numeric(obj_dim)-10^(10)),apply(sample_PF,2,max),c(1:length(sample_PF[,1])))
			#cell2<- create_cells(sample_PF,apply(sample_PF,2,min),apply(sample_PF,2,max),c(1:length(sample_PF[,1])))
			cell3<- create_cells(sample_PF[dist_id,],(numeric(obj_dim)-10^(10)),apply(sample_PF[dist_id,],2,max),c(1:length(sample_PF[dist_id,1])))
			
			



				for(id in 1:total_input_num)
				{
				#af1<-compute_PFES(RFM_post_mean_mat[id,],RFM_post_var_mat[id,],cell3)
				af1<-compute_PFES(post_mean_mat[id,],post_var_mat[id,],cell3)

				AF[id]<- AF[id]+af1#*(max(apply(t(sample_PF[dist_id,])-SAMPLE_Y_MAT[id,],2,min))>=0)
				#AF2[id]<- AF2[id]+compute_PFES2(RFM_post_mean_mat[id,],RFM_post_var_mat[id,],cell3)

				}

			}

		
			if(max(AF)==0)
			{
			new_X_num<- sample(c(1:total_input_num),1)
			} else {
			new_X_num<- order(AF,decreasing=T)[1]
			}

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}




		if(AF_num==11)
		{
		starttime<-Sys.time()

		distance_mat<- diag(ww)*0
			for(www in 1:ww)
			{
			distance_mat[www,]<- score_function(apply(t(Y)-Y[www,],2,min))
			}

		f.obj<-numeric(ww)+1
		f.con<-rbind(distance_mat,diag(ww))
		f.dir<-c()

			for(www in 1:(2*ww))
			{
			f.dir<-c(f.dir,">=")
			}

		f.rhs<- c((numeric(ww)+0.5),numeric(ww))
		
		lp_result<- lp("min",f.obj,f.con,f.dir,f.rhs)
		
		alpha_solution<- lp_result$solution

		gt_y<- numeric(total_input_num)
			
			for(www in 1:total_input_num)
			{
			gt_y[www] <- 1-2*sum(alpha_solution*score_function(apply(t(Y)-post_mean_mat[www,],2,min)))
			}


		deri_num_mat<-matrix(numeric(total_input_num*obj_dim*input_dim),total_input_num*obj_dim,input_dim)
		

			if(ww==1)
			{
			kvecd1<-matrix(numeric(total_input_num*obj_dim*iteration_num),total_input_num*obj_dim,iteration_num)
			kvecd2<-matrix(numeric(total_input_num*obj_dim*iteration_num),total_input_num*obj_dim,iteration_num)

			}


			for(www in 1:obj_dim)
			{	
			invC<-solve(cbind(Kvec_mat[((www-1)*total_input_num+X_num),1:ww])+diag(ww)*sigman)

			kvecd1[((www-1)*total_input_num+1):(www*total_input_num),ww]<- apply(cbind(all_input,t(new_X)%x%(numeric(total_input_num)+1),www),1,compute_kernel_deri1)
			kvecd2[((www-1)*total_input_num+1):(www*total_input_num),ww]<- apply(cbind(all_input,t(new_X)%x%(numeric(total_input_num)+1),www),1,compute_kernel_deri2)

			deri_num_mat[((www-1)*total_input_num+1):(www*total_input_num),]<-cbind(kvecd1[((www-1)*total_input_num+1):(www*total_input_num),1:ww]%*%invC%*%Y[,www],kvecd2[((www-1)*total_input_num+1):(www*total_input_num),1:ww]%*%invC%*%Y[,www])
			}

		deri_gt_mu<-matrix(numeric(total_input_num*input_dim),total_input_num,input_dim)

			for(www in 1:total_input_num)
			{
			dsf<-deri_score_function(apply(t(Y)-post_mean_mat[www,],2,min))
			number<-apply(t(Y)-post_mean_mat[www,],2,order)[1,]
			
			aabbcc<- deri_num_mat[((c(1:obj_dim)-1)*total_input_num+www),]
			aaabbbccc<-rbind(aabbcc[number,])

			deri_gt_mu[www,]<- -2*colSums(alpha_solution*dsf*aaabbbccc)

			}
		eta<-sqrt(rowSums(deri_gt_mu*deri_gt_mu))
		

		AF<-gt_y + beta*eta*sqrt(post_var_mat[,1])
			
	

		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==12)
		{
		starttime<-Sys.time()

		scalorization_vec<- rnorm(obj_dim,0,1) 
		scalorization_vec<-abs(scalorization_vec)/sum(abs(scalorization_vec))

		AF<- as.vector(ucb_mat%*%scalorization_vec)
		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}

		if(AF_num==13)
		{
		starttime<-Sys.time()

		scalorization_vec<- rnorm(obj_dim,0,1) 
		scalorization_vec<-abs(scalorization_vec)/sum(abs(scalorization_vec))

		sample_y_mat<-matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
			

			for(www in 1:obj_dim)
			{
			RFM_X<- RFM_mat[((www-1)*total_input_num+1):(www*total_input_num),]
			current_RFM_X<-rbind(RFM_X[X_num,])
			inv_tXX_I<- solve(t(current_RFM_X)%*%current_RFM_X+sigman*diag(RFM_sample_num))
			post_coefficient<-inv_tXX_I%*%t(current_RFM_X)%*%Y[,www]
			post_RFM_cov<- sigman*inv_tXX_I
			
			generated_coeficient<-t(rbind(mvrnorm(1,post_coefficient,(post_RFM_cov+10^(-8)*diag(RFM_sample_num)))))
			sample_y_mat[,www]<-RFM_X%*%generated_coeficient

			}



		AF<- as.vector(sample_y_mat%*%scalorization_vec)
		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==14)
		{
		starttime<-Sys.time()

		cPF_num<-compute_PF(Y)[,(obj_dim+1)]
		AF<-numeric(total_input_num)

		for(www in 1:total_input_num)
		{
		sample_mat<- mvrnorm(EMMDI_sample_num,post_mean_mat[www,],diag(post_var_mat[www,]))
		AF[www]<-compute_EMMDI(sample_mat,rbind(Y[cPF_num,]))
		}
		new_X_num<- order(AF,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}


	######################################################################

	new_X<- all_input[new_X_num,]
	new_Y<- true_function_mat[new_X_num,]+rnorm(obj_dim,0,sqrt(sigman))
		
	}

}


write.table(simple_PHV_regret_mat,"simple_PHV_regret_mat.txt")
write.table(simple_e_accuracy_mat,"simple_e_accuracy_mat.txt")
write.table(simple_e_accuracy2_mat,"simple_e_accuracy2_mat.txt")


write.table(inference_PHV_regret_mat,"inference_PHV_regret_mat.txt")
write.table(inference_e_accuracy_mat,"inference_e_accuracy_mat.txt")
write.table(inference_e_accuracy2_mat,"inference_e_accuracy2_mat.txt")


write.table(AF_computation_time_mat,"AF_computation_time_mat.txt")
write.table(X_num_mat,"X_num_mat.txt")




