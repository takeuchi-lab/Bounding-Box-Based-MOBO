library(nsga2R)
library(MASS)
library(lpSolve)

seed<-1
total_seed_num<-920
iteration_num<-500
setting_ID<- 3	#1:max, 2:min, 3:mean, 4:DR_mean
AF_num<-1

sigmaf1<-25
sigmaf2<-25
sigman<-10^(-6)
obj_dim<-2
beta<-3

EHVI_sample_num<-20

upper_diff_f<- 0.25

Y_mat<-read.csv("Y.csv")
Y_mat<-Y_mat[,2:3]
Y_mat<- as.matrix(Y_mat)

info<- read.csv("information.csv")
info<- as.matrix(info[,2:3])






num_diff_st<-numeric(429)

for(i in 1:429)
{
num_diff_st[i]<- sum(info[,1]==i)
}

id_num_diff_st<- cbind(c(1:429),num_diff_st)

for(i in 1:429)
{
info[info[,1]==i,2]<- c(1:num_diff_st[i])
}

modified_Y<- Y_mat-t(c(5,5))%x%rep(1,970)
true_Y_mat<- modified_Y[info[,2]<=10,]
modified_Y<-true_Y_mat

all_x<- as.matrix(read.csv("input_X.csv")[,-1])
all_x<- scale(all_x)

all_w<- as.matrix(read.csv("input_W.csv")[,-1])
all_w<- scale(all_w)

all_x_w<- cbind(all_x,all_w)

num_diff_st[num_diff_st>=10]<-10
all_x<-all_x[info[,2]<=10,]
all_x<- scale(all_x)
all_w<-all_w[info[,2]<=10,]
all_w<- scale(all_w)
all_x_w<- cbind(all_x,all_w)
info<- info[info[,2]<=10,]

id_num_diff_st<- cbind(c(1:429),num_diff_st)

total_input_num<-length(all_x[,1])



if(setting_ID==1)
{
max_Y_mat<- matrix(numeric(429*2),429,2)
for(i in 1:429)
{
max_Y_mat[i,]<-  apply(rbind(true_Y_mat[info[,1]==i,]),2,max)


}
true_risk_mat<- max_Y_mat
}

if(setting_ID==2)
{
min_Y_mat<- matrix(numeric(429*2),429,2)
for(i in 1:429)
{
min_Y_mat[i,]<-  apply(rbind(true_Y_mat[info[,1]==i,]),2,min)


}
true_risk_mat<- min_Y_mat

}


if(setting_ID==3)
{
mean_Y_mat<- matrix(numeric(429*2),429,2)
for(i in 1:429)
{
mean_Y_mat[i,]<-  apply(rbind(true_Y_mat[info[,1]==i,]),2,mean)


}
true_risk_mat<- mean_Y_mat

}




compute_DR<-function(x)
{
	if(length(x)==1)
	{
	result<-x
	} else {
	ref_p<-  rep(1,length(x))/length(x)
	f_ref_p<-sum(x*ref_p)
	grid_num_w<- length(x)


		f1<- (numeric(grid_num_w)+1)%x%c(1,-1)
		D<-diag(2*grid_num_w)
		f2<-numeric(2*grid_num_w)+1
		f3<- diag(grid_num_w)%x% t(c(1,-1))

		f.con <- rbind(f1,D,f2,f3,f3)


		f.dir <- c("==")

		for(i in 1:(2*grid_num_w) )
		{
		f.dir <- c(f.dir,">=")
		}

		f.dir <- c(f.dir,"<=")

		for(i in 1:grid_num_w )
		{
		f.dir <- c(f.dir,"<=")
		}

		for(i in 1:grid_num_w )
		{
		f.dir <- c(f.dir,">=")
		}


		f.rhs <- c(numeric((2*grid_num_w)+1),upper_diff_f,(numeric(grid_num_w)+1-ref_p),(numeric(grid_num_w)-ref_p))
		f.rhs_c <- c(numeric((2*grid_num_w)+1),upper_diff_f,(numeric(grid_num_w)+1-ref_p),(numeric(grid_num_w)-ref_p))

		f.obj <-x%x%c(1,-1)


resu<-lp ("min", f.obj, f.con, f.dir, f.rhs)
kekka<-resu$objval


result<-kekka+f_ref_p

	}
result
}









if(setting_ID==4)
{
DR_mean_Y_mat<- matrix(numeric(429*2),429,2)




for(i in 1:429)
{
DR_mean_Y_mat[i,]<-  apply(rbind(true_Y_mat[info[,1]==i,]),2,compute_DR)


}
true_risk_mat<- DR_mean_Y_mat

}

#########################################

all_x_w2<-c()
for(i in 1:429)
{
all_x_w2<-rbind(all_x_w2,t(colMeans(rbind(all_x_w[info[,1]==i,])))%x%rep(1,sum(info[,1]==i)))
}

all_x<- all_x_w2
all_w<- all_x_w

all_x<- scale(all_x)
all_w<- scale(all_w)

all_x_w<- cbind(all_x,all_w)

x_k_mat<- matrix(numeric(total_input_num*total_input_num),total_input_num,total_input_num)
w_k_mat<- matrix(numeric(total_input_num*total_input_num),total_input_num,total_input_num)
x_w_k_mat<- matrix(numeric(total_input_num*total_input_num),total_input_num,total_input_num)

for(i in 1:total_input_num)
{
	for(j in 1:total_input_num)
	{
	x_k_mat[i,j]<- sum((all_x[i,]-all_x[j,])^2)
	w_k_mat[i,j]<- sum((all_w[i,]-all_w[j,])^2)
	x_w_k_mat[i,j]<- sum((all_x_w[i,]-all_x_w[j,])^2)

	}
}

Len1<-c()
Len2<-c()
Len3<-c()

for(i in 1:(total_input_num-1))
{
	Len1<-c(Len1,x_k_mat[i,(i+1):total_input_num])
	Len2<-c(Len2,w_k_mat[i,(i+1):total_input_num])
	Len3<-c(Len3,x_w_k_mat[i,(i+1):total_input_num])

}

length1<- median(Len1)
length2<- median(Len2)
length3<- median(Len3)


sigmaf<- c(sigmaf1,sigmaf2)

compute_kernel<- function(a)
{
x1<- a[1:51]
w1<- a[52:102]
x2<- a[103:153]
w2<- a[154:204]
l<-a[205]

d1<- sum((x1-x2)^2)
d2<- sum((w1-w2)^2)

sigmaf[l]*exp(-d1/length_scale[1])*exp(-d2/length_scale[2])

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


standard_normal_pdf<-function(x)
{
(2*pi)^(-1/2)*exp(-x^2/2)
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



###############################################################
##################     Pareto Front     #######################
###############################################################

true_PF_PFn<-compute_PF(true_risk_mat)
true_PF<-true_PF_PFn[,1:obj_dim]
true_PFn<-true_PF_PFn[,(obj_dim+1)]

ref_point<- apply(true_risk_mat,2,min)
max_point<- apply(true_risk_mat,2,max)


#set.seed(1)
true_PF_cell<- create_cells(true_PF,ref_point,max_point,c(1:length(true_PF[,1])))
true_PHV<- sum(apply(rbind(rbind(true_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(true_PF_cell)[,1:(obj_dim)]),1,prod))

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
Kvec_mat_tilde<-matrix(numeric(429*iteration_num),429,iteration_num)%x%numeric(obj_dim)


post_mean_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
post_var_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)

post_mean_mat_tilde<- matrix(numeric(429*obj_dim),429,obj_dim)
post_var_mat_tilde<- matrix(numeric(429*obj_dim),429,obj_dim)

ucb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
lcb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
post_cov_mat<-diag(total_input_num)%x%numeric(obj_dim)

svd_mat<-diag(total_input_num)%x%numeric(obj_dim)

#original_all_x_w<-all_x_w

#all_x_w<-all_x_w*10

#length_scale<-c(length1,length2)/2
length_scale<-c(length3,length3)/2

prior_var_tilde1<-numeric(429)
prior_var_tilde2<-numeric(429)

for(i in 1:429)
{

	ABC<-  cbind(rbind(all_x[info[,1]==i,]),rbind(all_w[info[,1]==i,]))
	ABC2<- cbind(ABC%x%rep(1,length(ABC[,1])),rep(1,length(ABC[,1]))%x%ABC)
	prior_var_tilde1[i]<-mean(apply(cbind(ABC2,1),1,compute_kernel))
	prior_var_tilde2[i]<-mean(apply(cbind(ABC2,2),1,compute_kernel))

}

prior_var_tilde<-cbind(prior_var_tilde1,prior_var_tilde2)
all_computation_time_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)


compute_tilde<-function(a)
{
n<-c()
	for(i in 1:429)
	{
		n[i]<-mean(a[info[,1]==i])
	}
n
}

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


new_X_num<- (seed-1)*total_seed_num+w
new_X<- all_x_w[new_X_num,]
new_Y<- true_Y_mat[new_X_num,]#+rnorm(obj_dim,0,sqrt(sigman))
	

	for(ww in 1:iteration_num)
	{
all_comp_start_time<-Sys.time()

	X<-rbind(X,new_X)	
	X_num<-c(X_num,new_X_num)
	Y<-rbind(Y,new_Y)
	G1G2<- diag(ww)%x%rep(1,2)
		for(www in 1:obj_dim)
		{
		Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),ww]<-t(matrix(apply(cbind(all_x_w,t(new_X)%x%(numeric(total_input_num)+1),www),1,compute_kernel),1,total_input_num))

		Kvec_mat_tilde[((www-1)*429+1):(www*429),ww]<-compute_tilde(Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),ww])


		partial_Kvec_mat<-cbind(Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),1:ww])
		partial_Kvec_mat_tilde<-cbind(Kvec_mat_tilde[((www-1)*429+1):(www*429),1:ww])
			if(www==1)
			{
			G1G2[1:ww,]<-partial_Kvec_mat[X_num,]
			} else {
			G1G2[(ww+1):(2*ww),]<-partial_Kvec_mat[X_num,]
			}
		Cinv<-solve(partial_Kvec_mat[X_num,]+sigman*diag(ww))


		

		post_mean_mat[,www]<-as.vector(partial_Kvec_mat%*%Cinv%*%Y[,www])
		post_var_mat[,www]<- abs(sigmaf[www]-rowSums((partial_Kvec_mat%*%Cinv)*partial_Kvec_mat))

		post_mean_mat_tilde[,www]<-as.vector(partial_Kvec_mat_tilde%*%Cinv%*%Y[,www])
		post_var_mat_tilde[,www]<- abs(prior_var_tilde[,www]-rowSums((partial_Kvec_mat_tilde%*%Cinv)*partial_Kvec_mat_tilde))



		}
	G1<- G1G2[1:ww,]
	G2<- G1G2[(ww+1):(2*ww),]

	ucb_mat<- post_mean_mat+beta*sqrt(post_var_mat)
	lcb_mat<- post_mean_mat-beta*sqrt(post_var_mat)

	predicted_PF_PFn<-compute_PF(rbind(true_risk_mat[info[X_num,1],]))
	predicted_PF<-rbind(predicted_PF_PFn[,1:obj_dim])
	predicted_PFn<-X_num[predicted_PF_PFn[,(obj_dim+1)]]

predicted_PF_cell<-create_cells(predicted_PF,ref_point,apply(predicted_PF,2,max),c(1:length(predicted_PF[,1])))

	predicted_PF_PHV<-  sum(apply(rbind(predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-predicted_PF_cell[,(1):(obj_dim)]),1,prod))
	
	
	simple_PHV_regret<- c(simple_PHV_regret,max(true_PHV-predicted_PF_PHV,0))
	simple_PHV_regret_mat[w,ww]<-	min(simple_PHV_regret)

	simple_e_accuracy<- c(simple_e_accuracy,compute_e_accuracy(true_PF,predicted_PF))
	simple_e_accuracy_mat[w,ww]<-	min(simple_e_accuracy)

	simple_e_accuracy2<- c(simple_e_accuracy2,compute_e_accuracy2(true_PF,predicted_PF))
	simple_e_accuracy2_mat[w,ww]<-	min(simple_e_accuracy2)

	


	#if(setting_ID==3)
	#{
	#	compute_LCB<-function(i)
	#	{
	#	apply(rbind(lcb_mat[info[,1]==i,]),2,mean)
	#	}
	#
	#	compute_UCB<-function(i)
	#	{
	#	apply(rbind(ucb_mat[info[,1]==i,]),2,mean)
	#	}
	#}

	




	#LCB_mat<-t(apply(cbind(c(1:429)),1,compute_LCB))
	#UCB_mat<-t(apply(cbind(c(1:429)),1,compute_UCB))

	pPF2<-rbind(compute_PF(rbind(post_mean_mat_tilde[info[X_num,1],]))[,1:(obj_dim)])
	pPF2_num<-info[X_num[compute_PF(rbind(post_mean_mat_tilde[info[X_num,1],]))[,(obj_dim+1)]],1]
	reported_Y<-rbind(true_risk_mat[pPF2_num,])
	inference_predicted_PF_PFn<-compute_PF(rbind(true_risk_mat[pPF2_num,]))
	inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])
	inference_predicted_PFn<-pPF2_num[inference_predicted_PF_PFn[,(obj_dim+1)]]

	inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
	inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))

	pPF2_PF_cell<- create_cells(pPF2,ref_point,apply(pPF2,2,max),c(1:length(pPF2[,1])))
	pPF2_PHV<-  sum(apply(rbind(pPF2_PF_cell[,(obj_dim+1):(2*obj_dim)]-pPF2_PF_cell[,(1):(obj_dim)]),1,prod))



	inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
	inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret[ww]

	inference_e_accuracy<- c(inference_e_accuracy,compute_e_accuracy(true_PF,inference_predicted_PF))
	inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy[ww]

	inference_e_accuracy2<- c(inference_e_accuracy2,compute_e_accuracy2(true_PF,reported_Y))
	inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2[ww]

	X_num_mat[w,ww]<-new_X_num

	##########################################
		if(AF_num==1)
		{
		starttime<-Sys.time()

		compute_EI<-function(x)
		{
		mu<- post_mean_mat_tilde[x,]
		var<- post_var_mat_tilde[x,]
		sample<-  mvrnorm(EHVI_sample_num,mu,diag(var))
		result<- numeric(EHVI_sample_num)

			for(i in 1:EHVI_sample_num)
			{
			c<-compute_e_accuracy(rbind(sample[i,]),pPF2)
				if(c<=0)
				{
				result[i]<-0
				} else {
				D<- rbind(pPF2,sample[i,])
				D_cell<-create_cells(D,ref_point,apply(D,2,max),c(1:length(D[,1])))

				D_PF_PHV<-  sum(apply(rbind(D_cell[,(obj_dim+1):(2*obj_dim)]-D_cell[,(1):(obj_dim)]),1,prod))
				result[i]<-D_PF_PHV-pPF2_PHV
	

				}
			}
		mean(result)
		}

		AFX<- apply(cbind(c(1:429)),1,compute_EI)
		next_design<- order(AFX,decreasing=T)[1]
		new_input_design<-  rbind(all_x[info[,1]==next_design,])[1,]


		new_var_tilde_function<-function(iii)
			{
			W<- rbind(all_w[info[,1]==next_design,])[iii,]
			A<-rbind(X,c(new_input_design,W))
			B1<-apply(cbind(A,t(new_input_design)%x%rep(1,length(A[,1])),t(W)%x%rep(1,length(A[,1])),1),1,compute_kernel)
			B2<-apply(cbind(A,t(new_input_design)%x%rep(1,length(A[,1])),t(W)%x%rep(1,length(A[,1])),2),1,compute_kernel)

			new_G1<-   rbind(cbind(G1,B1[1:ww]),B1)
			new_G2<-   rbind(cbind(G2,B2[1:ww]),B2)
			
			new_C1<- new_G1+sigman*diag((ww+1))
			new_C2<- new_G2+sigman*diag((ww+1))

			new_C1_inv<- solve(new_C1)
			new_C2_inv<- solve(new_C2)


						
			W3<- rbind(all_w[info[,1]==next_design,])
			H1<-colMeans(matrix(apply(cbind(rep(1,length(A[,1]))%x%cbind(rbind(new_input_design)%x%rep(1,length(W3[,1])),W3) ,A%x%rep(1,length(W3[,1])),1),1,compute_kernel),length(W3[,1]),length(A[,1])))
			H2<-colMeans(matrix(apply(cbind(rep(1,length(A[,1]))%x%cbind(rbind(new_input_design)%x%rep(1,length(W3[,1])),W3) ,A%x%rep(1,length(W3[,1])),2),1,compute_kernel),length(W3[,1]),length(A[,1])))


			h1<- as.vector(t(H1)%*%new_C1_inv%*%H1)
			h2<- as.vector(t(H2)%*%new_C2_inv%*%H2)
			
			v<-abs(post_var_mat_tilde[next_design,]-c(h1,h2))
			sum(log(v))
			}

		AFW<- apply(cbind(c(1:sum(info[,1]==next_design))),1,new_var_tilde_function)
		next_env<- order(AFW)[1]

		new_X_num<-c(1:920)[(info[,1]==next_design)&(info[,2]==next_env)]



		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		}



	new_X<- all_x_w[new_X_num,]
	new_Y<- true_Y_mat[new_X_num,]#+rnorm(obj_dim,0,sqrt(sigman))
		all_comp_end_time<-Sys.time()
all_computation_time_mat[w,ww]<- difftime(all_comp_end_time, all_comp_start_time, units = "sec") 


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
write.table(all_computation_time_mat,"all_computation_time_mat.txt")


