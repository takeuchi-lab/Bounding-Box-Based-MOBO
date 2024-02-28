library(nsga2R)
library(MASS)
library(lpSolve)

seed<-1
total_seed_num<-920
iteration_num<-500
setting_ID<- 2	#1:max, 2:min, 3:mean, 4:DR_mean
AF_num<-3

sigmaf1<-25
sigmaf2<-25
sigman<-10^(-6)
obj_dim<-2
beta<-3

simulator_setting<-1

SABBa_AF_sample_num<-20
SABBa_new_input_sample<-10
SABBa_initial<-0.5#0.5
SABBa_ratio<-0.9
SABBa_min_ratio<-0.01


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

compute_kernel_SABBa<- function(a)
{
x1<- a[1:51]
x2<- a[52:102]

d1<- sum((x1-x2)^2)


sigmaf[1]*exp(-d1/length_scale_SABBa)

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
post_mean_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
post_var_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
ucb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
lcb_mat<- matrix(numeric(total_input_num*obj_dim),total_input_num,obj_dim)
post_cov_mat<-diag(total_input_num)%x%numeric(obj_dim)

svd_mat<-diag(total_input_num)%x%numeric(obj_dim)

#original_all_x_w<-all_x_w

#all_x_w<-all_x_w*10

length_scale<-c(length3,length3)/2

length_scale_SABBa<- length1/2

threshold_s1<-max(true_risk_mat[,1])-min(true_risk_mat[,1])
threshold_s2<-max(true_risk_mat[,2])-min(true_risk_mat[,2])

all_x_part<- all_x[,1:51]
all_x_part<- all_x_part[info[,2]==1,]

id_num_diff_st2<- cbind(id_num_diff_st,id_num_diff_st)
for(i in 1:429)
{
id_num_diff_st2[i,3]<- min(c(1:920)[info[,1]==i])
id_num_diff_st2[i,4]<- max(c(1:920)[info[,1]==i])

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

all_computation_time_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)



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
	


	X_tilde_id<-sample(c(1:429),SABBa_new_input_sample)
	X_tilde_id<- sort(X_tilde_id)
	ratio<-max(SABBa_initial,SABBa_min_ratio)
	th1<- threshold_s1*ratio
	th2<- threshold_s2*ratio
	new_X_tilde_generate_counter<-0

	outer_GP_mu<-c()
	outer_GP_var<-c()
	outer_GP_id<-c()


	for(ww in 1:iteration_num)
	{
all_comp_start_time<-Sys.time()

	X<-rbind(X,new_X)	
	X_num<-c(X_num,new_X_num)
	Y<-rbind(Y,new_Y)

		for(www in 1:obj_dim)
		{
		Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),ww]<-t(matrix(apply(cbind(all_x_w,t(new_X)%x%(numeric(total_input_num)+1),www),1,compute_kernel),1,total_input_num))

		partial_Kvec_mat<-cbind(Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),1:ww])
		Cinv<-solve(partial_Kvec_mat[X_num,]+1000*sigman*diag(ww))

		post_mean_mat[,www]<-as.vector(partial_Kvec_mat%*%Cinv%*%Y[,www])
		post_var_mat[,www]<- abs(sigmaf[www]-rowSums((partial_Kvec_mat%*%Cinv)*partial_Kvec_mat))

		}

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

	








##############################################################################################
		if(setting_ID==2)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)

			if(new_X_tilde_generate_counter==1)
			{
			ratio<-max(ratio*SABBa_ratio,SABBa_min_ratio)
			th1<-threshold_s1*ratio
			th2<-threshold_s2*ratio
			X_tilde_id<-c(X_tilde_id,sample(c(1:429),SABBa_new_input_sample))
			X_tilde_id<- unique(X_tilde_id)
			X_tilde_id<- sort(X_tilde_id)
			new_X_tilde_generate_counter<-0
			}


		epsilon_mat<- ucb_mat-post_mean_mat
		risk_epsilon_mat<- matrix(numeric(429*2),429,2)
		risk_post_mean_mat<-matrix(numeric(429*2),429,2)

			for(i in 1:429)
			{
			risk_post_mean_mat[i,]<-apply(rbind(post_mean_mat[id_num_diff_st2[i,3]:id_num_diff_st2[i,4],]),2,min)
			risk_epsilon_mat[i,]<-apply(rbind(abs(apply(rbind(ucb_mat[id_num_diff_st2[i,3]:id_num_diff_st2[i,4],]),2,min)-risk_post_mean_mat[i,]),abs(apply(rbind(lcb_mat[id_num_diff_st2[i,3]:id_num_diff_st2[i,4],]),2,min)-risk_post_mean_mat[i,])),2,max)


			}

		X_tilde_risk_post_mean_mat<- risk_post_mean_mat[X_tilde_id,]
		X_tilde_risk_epsilon_mat<- risk_epsilon_mat[X_tilde_id,]
		X_tilde_risk_LCB_mat<-X_tilde_risk_post_mean_mat-X_tilde_risk_epsilon_mat
		X_tilde_risk_UCB_mat<-X_tilde_risk_post_mean_mat+X_tilde_risk_epsilon_mat

			if( sum((X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2))>0)
			{
			outer_GP_id<-c(outer_GP_id,X_tilde_id[(X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2)])
			outer_GP_mu<- rbind(outer_GP_mu,rbind(X_tilde_risk_post_mean_mat[(X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2),]))
			outer_GP_var<- rbind(outer_GP_var,(rbind(X_tilde_risk_epsilon_mat[(X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2),])/3)^2)

				sort_oGPid<-unique(outer_GP_id)
				nsort_oGPid<-sort_oGPid
				for(aabb in 1:length(sort_oGPid))
				{
				nsort_oGPid[aabb]<- max(c(1:length(outer_GP_id))[outer_GP_id==sort_oGPid[aabb]])
				}
			outer_GP_id<- outer_GP_id[nsort_oGPid]
			outer_GP_mu<-rbind(outer_GP_mu[nsort_oGPid,])
			outer_GP_var<-rbind(outer_GP_var[nsort_oGPid,])
			}

			if(length(outer_GP_id)>0)
			{
			outer_gram_mat<- matrix(apply(cbind(rbind(all_x_part[outer_GP_id,])%x%rep(1,length(outer_GP_id)),rep(1,length(outer_GP_id))%x%rbind(all_x_part[outer_GP_id,])),1,compute_kernel_SABBa),length(outer_GP_id),length(outer_GP_id))
				if(length(outer_GP_id)==1)
				{
				outer_GP_C1<- outer_gram_mat+outer_GP_var[1,1]+10^(-3)
				outer_GP_C2<- outer_gram_mat+outer_GP_var[1,2]+10^(-3)

				} else {
				outer_GP_C1<- outer_gram_mat+diag(outer_GP_var[,1]+10^(-3))
				outer_GP_C2<- outer_gram_mat+diag(outer_GP_var[,2]+10^(-3))
				}
			outer_GP_C1_inv<- solve(outer_GP_C1)
			outer_GP_C2_inv<- solve(outer_GP_C2)

			outer_GP_kvec<- t(matrix(apply(cbind(all_x_part%x%rep(1,length(outer_GP_id)),rep(1,429)%x%rbind(all_x_part[outer_GP_id,])),1,compute_kernel_SABBa),length(outer_GP_id),429))
			outer_GP_risk_post_mean_mat<- cbind(outer_GP_kvec%*%outer_GP_C1_inv%*%outer_GP_mu[,1],outer_GP_kvec%*%outer_GP_C2_inv%*%outer_GP_mu[,2])
			outer_GP_risk_post_var_mat<- cbind(abs(sigmaf[1]-rowSums((outer_GP_kvec%*%outer_GP_C1_inv)*outer_GP_kvec)),abs(sigmaf[1]-rowSums((outer_GP_kvec%*%outer_GP_C2_inv)*outer_GP_kvec)))
			outer_GP_risk_epsilon_mat<- 3*sqrt(outer_GP_risk_post_var_mat)
			outer_GP_X_tilde_risk_epsilon_mat<-rbind(outer_GP_risk_epsilon_mat[X_tilde_id,])
			}

		BB_box_dominate<- compute_dist(X_tilde_risk_UCB_mat,	rbind(compute_PF_2dim(X_tilde_risk_LCB_mat)[,1:2]))

		check1<- (BB_box_dominate<0)
		check2<- (X_tilde_risk_epsilon_mat[,1]<th1)&(X_tilde_risk_epsilon_mat[,2]<th2)
			if(length(outer_GP_id)>0)
			{
			check3<- (outer_GP_X_tilde_risk_epsilon_mat[,1]<th1)&(outer_GP_X_tilde_risk_epsilon_mat[,2]<th2)
			}
		
			if(length(outer_GP_id)>0)
			{
			check<- (1-check1)*(1-check2)*(1-check3)
			} else {
			check<- (1-check1)*(1-check2)
			}



			if(sum(check)==0)
			{

			new_design_num<-sample(X_tilde_id,1)

			new_X_tilde_generate_counter<-1

			} else {
				if(sum(check)==1)
				{
				new_design_num<-X_tilde_id[check==1]

				} else {
				B_minus<- X_tilde_risk_LCB_mat
				B_plus<- X_tilde_risk_UCB_mat




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
			AF<- apply( cbind(c(1:length(X_tilde_id))),1,function_pop)
			AF<-AF*check
			new_design_num<-X_tilde_id[order(AF,decreasing=T)[1]]


				}
			}

			if(id_num_diff_st2[new_design_num,2]==1)
			{
			new_env_num<-1
			} else {

			check_d<-c(1:920)[info[,1]==new_design_num]
			check_d2<-check_d*0
				for(i in 1:length(check_d2))
				{
				check_d2[i]<- sum(X_num==check_d[i])
				}

				if(sum(check_d2)==0)
				{
				current_m<-apply(post_mean_mat[id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4],],2,max)

				} else {
					check_3<-check_d[check_d2>0]
					c4<-matrix(numeric(length(check_3)*2),length(check_3),2)
					for(i in 1:length(check_3))
					{
					c4[i,]<-rbind(Y[X_num==check_3[i],])[1,]
					}
				current_m<- apply(c4,2,min)
				}
			c_k1<- (current_m[1]-lcb_mat[id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4],1])/(2*epsilon_mat[id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4],1])
			c_k1<- c_k1*(c_k1>0)

			c_k2<- (current_m[2]-lcb_mat[id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4],2])/(2*epsilon_mat[id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4],2])
			c_k2<- c_k2*(c_k2>0)



			c_bar_1<- (c_k1-min(c_k1))/(max(c_k1)-min(c_k1))
			c_bar_2<- (c_k2-min(c_k2))/(max(c_k2)-min(c_k2))
			wei<-(risk_epsilon_mat[new_design_num,]/th2)^2
			AFW<- (wei[1]*c_bar_1+wei[2]*c_bar_2)
			c6<-id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4]
			for(i in 1:length(c6))
			{
			c6[i]<- sum(X_num==c6[i])
			}

			 AFW[c6>0]<-  -Inf

			new_env_num<- order(AFW,decreasing=T)[1]

			}

		


		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 


			if(ww==1)
			{
			inference_predicted_PF_PFn<-compute_PF_2dim(X_tilde_risk_post_mean_mat)[,3]

			reported_Y_mat<-rbind(true_risk_mat[X_tilde_id[inference_predicted_PF_PFn],])

			inference_predicted_PF_PFn<-compute_PF_2dim(reported_Y_mat)
			inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


			inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
			inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


			#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
			inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

			#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
			inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

			#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
			inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y_mat),0)
			} 



			if((new_X_tilde_generate_counter==0)&(ww!=1))
			{
			inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret_mat[w,(ww-1)]

			inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy_mat[w,(ww-1)]

			inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2_mat[w,(ww-1)]


			} else {
			inference_predicted_PF_PFn<-compute_PF_2dim(X_tilde_risk_post_mean_mat)[,3]

			reported_Y_mat<-rbind(true_risk_mat[X_tilde_id[inference_predicted_PF_PFn],])

			inference_predicted_PF_PFn<-compute_PF_2dim(reported_Y_mat)
			inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


			inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
			inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


			#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
			inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

			#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
			inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

			#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
			inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y_mat),0)

			}

			

		
		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		


		

		}
#############################################################################################################

		if(setting_ID==3)
		{
		starttime<-Sys.time()

		set.seed(((seed-1)*total_seed_num+w)*1000+ww)

			if(new_X_tilde_generate_counter==1)
			{
			ratio<-max(ratio*SABBa_ratio,SABBa_min_ratio)
			th1<-threshold_s1*ratio
			th2<-threshold_s2*ratio
			X_tilde_id<-c(X_tilde_id,sample(c(1:429),SABBa_new_input_sample))
			X_tilde_id<- unique(X_tilde_id)
			X_tilde_id<- sort(X_tilde_id)
			new_X_tilde_generate_counter<-0
			}


		epsilon_mat<- ucb_mat-post_mean_mat
		risk_epsilon_mat<- matrix(numeric(429*2),429,2)
		risk_post_mean_mat<-matrix(numeric(429*2),429,2)

			for(i in 1:429)
			{
			risk_epsilon_mat[i,]<-colMeans(rbind(epsilon_mat[id_num_diff_st2[i,3]:id_num_diff_st2[i,4],]))
			risk_post_mean_mat[i,]<-colMeans(rbind(post_mean_mat[id_num_diff_st2[i,3]:id_num_diff_st2[i,4],]))

			}

		X_tilde_risk_post_mean_mat<- risk_post_mean_mat[X_tilde_id,]
		X_tilde_risk_epsilon_mat<- risk_epsilon_mat[X_tilde_id,]
		X_tilde_risk_LCB_mat<-X_tilde_risk_post_mean_mat-X_tilde_risk_epsilon_mat
		X_tilde_risk_UCB_mat<-X_tilde_risk_post_mean_mat+X_tilde_risk_epsilon_mat

			if( sum((X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2))>0)
			{
			outer_GP_id<-c(outer_GP_id,X_tilde_id[(X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2)])
			outer_GP_mu<- rbind(outer_GP_mu,rbind(X_tilde_risk_post_mean_mat[(X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2),]))
			outer_GP_var<- rbind(outer_GP_var,(rbind(X_tilde_risk_epsilon_mat[(X_tilde_risk_epsilon_mat[,1]<th1)& (X_tilde_risk_epsilon_mat[,2]<th2),])/3)^2)

				sort_oGPid<-unique(outer_GP_id)
				nsort_oGPid<-sort_oGPid
				for(aabb in 1:length(sort_oGPid))
				{
				nsort_oGPid[aabb]<- max(c(1:length(outer_GP_id))[outer_GP_id==sort_oGPid[aabb]])
				}
			outer_GP_id<- outer_GP_id[nsort_oGPid]
			outer_GP_mu<-rbind(outer_GP_mu[nsort_oGPid,])
			outer_GP_var<-rbind(outer_GP_var[nsort_oGPid,])
			}

			if(length(outer_GP_id)>0)
			{
			outer_gram_mat<- matrix(apply(cbind(rbind(all_x_part[outer_GP_id,])%x%rep(1,length(outer_GP_id)),rep(1,length(outer_GP_id))%x%rbind(all_x_part[outer_GP_id,])),1,compute_kernel_SABBa),length(outer_GP_id),length(outer_GP_id))
				if(length(outer_GP_id)==1)
				{
				outer_GP_C1<- outer_gram_mat+outer_GP_var[1,1]+10^(-6)
				outer_GP_C2<- outer_gram_mat+outer_GP_var[1,2]+10^(-6)

				} else {
				outer_GP_C1<- outer_gram_mat+diag(outer_GP_var[,1]+10^(-6))
				outer_GP_C2<- outer_gram_mat+diag(outer_GP_var[,2]+10^(-6))
				}
			outer_GP_C1_inv<- solve(outer_GP_C1)
			outer_GP_C2_inv<- solve(outer_GP_C2)

			outer_GP_kvec<- t(matrix(apply(cbind(all_x_part%x%rep(1,length(outer_GP_id)),rep(1,429)%x%rbind(all_x_part[outer_GP_id,])),1,compute_kernel_SABBa),length(outer_GP_id),429))
			outer_GP_risk_post_mean_mat<- cbind(outer_GP_kvec%*%outer_GP_C1_inv%*%outer_GP_mu[,1],outer_GP_kvec%*%outer_GP_C2_inv%*%outer_GP_mu[,2])
			outer_GP_risk_post_var_mat<- cbind(abs(sigmaf[1]-rowSums((outer_GP_kvec%*%outer_GP_C1_inv)*outer_GP_kvec)),abs(sigmaf[1]-rowSums((outer_GP_kvec%*%outer_GP_C2_inv)*outer_GP_kvec)))
			outer_GP_risk_epsilon_mat<- 3*sqrt(outer_GP_risk_post_var_mat)
			outer_GP_X_tilde_risk_epsilon_mat<-rbind(outer_GP_risk_epsilon_mat[X_tilde_id,])
			}

		BB_box_dominate<- compute_dist(X_tilde_risk_UCB_mat,	rbind(compute_PF_2dim(X_tilde_risk_LCB_mat)[,1:2]))

		check1<- (BB_box_dominate<0)
		check2<- (X_tilde_risk_epsilon_mat[,1]<th1)&(X_tilde_risk_epsilon_mat[,2]<th2)
			if(length(outer_GP_id)>0)
			{
			check3<- (outer_GP_X_tilde_risk_epsilon_mat[,1]<th1)&(outer_GP_X_tilde_risk_epsilon_mat[,2]<th2)
			}
		
			if(length(outer_GP_id)>0)
			{
			check<- (1-check1)*(1-check2)*(1-check3)
			} else {
			check<- (1-check1)*(1-check2)
			}



			if(sum(check)==0)
			{

			new_design_num<-sample(X_tilde_id,1)

			new_X_tilde_generate_counter<-1

			} else {
				if(sum(check)==1)
				{
				new_design_num<-X_tilde_id[check==1]

				} else {
				B_minus<- X_tilde_risk_LCB_mat
				B_plus<- X_tilde_risk_UCB_mat




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
			AF<- apply( cbind(c(1:length(X_tilde_id))),1,function_pop)
			AF<-AF*check
			new_design_num<-X_tilde_id[order(AF,decreasing=T)[1]]


				}
			}

			if(id_num_diff_st2[new_design_num,2]==1)
			{
			new_env_num<-1
			} else {

			A1<-epsilon_mat[id_num_diff_st2[new_design_num,3]:id_num_diff_st2[new_design_num,4],]
			A2<-apply(A1,2,max)-apply(A1,2,min)
			c_bar_1<- (A1[,1]-min(A1[,1]))/A2[1]
			c_bar_2<- (A1[,2]-min(A1[,2]))/A2[2]
			wei<-(risk_epsilon_mat[new_design_num,]/th2)^2
			AFW<- (wei[1]*c_bar_1+wei[2]*c_bar_2)
			new_env_num<- order(AFW,decreasing=T)[1]

			}

		


		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 


			if(ww==1)
			{
			inference_predicted_PF_PFn<-compute_PF_2dim(X_tilde_risk_post_mean_mat)[,3]

			reported_Y_mat<-rbind(true_risk_mat[X_tilde_id[inference_predicted_PF_PFn],])

			inference_predicted_PF_PFn<-compute_PF_2dim(reported_Y_mat)
			inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


			inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
			inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


			#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
			inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

			#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
			inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

			#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
			inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y_mat),0)
			} 



			if((new_X_tilde_generate_counter==0)&(ww!=1))
			{
			inference_PHV_regret_mat[w,ww]<-	inference_PHV_regret_mat[w,(ww-1)]

			inference_e_accuracy_mat[w,ww]<-	inference_e_accuracy_mat[w,(ww-1)]

			inference_e_accuracy2_mat[w,ww]<-	inference_e_accuracy2_mat[w,(ww-1)]


			} else {
			inference_predicted_PF_PFn<-compute_PF_2dim(X_tilde_risk_post_mean_mat)[,3]

			reported_Y_mat<-rbind(true_risk_mat[X_tilde_id[inference_predicted_PF_PFn],])

			inference_predicted_PF_PFn<-compute_PF_2dim(reported_Y_mat)
			inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])


			inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
			inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


			#inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
			inference_PHV_regret_mat[w,ww]<-	max(true_PHV-inference_predicted_PF_PHV,0)

			#inference_e_accuracy<- c(inference_e_accuracy,max(compute_e_accuracy(true_PF,inference_predicted_PF),0))
			inference_e_accuracy_mat[w,ww]<-	max(compute_e_accuracy(true_PF,inference_predicted_PF),0)

			#inference_e_accuracy2<- c(inference_e_accuracy2,max(compute_e_accuracy2(true_PF,inference_predicted_PF),0))
			inference_e_accuracy2_mat[w,ww]<-	max(compute_e_accuracy2(true_PF,reported_Y_mat),0)

			}

			

		
		endtime<-Sys.time()

		AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 
		


		

		}
#############################################################################################################



		

		




	new_X_num<-c(1:920)[info[,1]==new_design_num & info[,2]==new_env_num][1]

	X_num_mat[w,ww]<-new_X_num



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


