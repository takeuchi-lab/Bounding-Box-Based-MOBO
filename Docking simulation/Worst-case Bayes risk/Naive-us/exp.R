library(nsga2R)
library(MASS)
library(lpSolve)

seed<-1
total_seed_num<-429
iteration_num<-500
setting_ID<- 4	#1:max, 2:min, 3:mean, 4:DR_mean
AF_num<-2	#1:random, 2:US, 3:EMmI, 4:ePAL

sigmaf1<-25
sigmaf2<-25
sigman<-10^(-6)
obj_dim<-2
beta<-3

total_input_num<-429
EMMDI_sample_num<-20

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



num_diff_st[num_diff_st>=10]<-10
all_x<-all_x[info[,2]<=10,]
all_x<- scale(all_x)
all_w<-all_w[info[,2]<=10,]
all_w<- scale(all_w)
all_x_w<- cbind(all_x,all_w)

info<- info[info[,2]<=10,]

id_num_diff_st<- cbind(c(1:429),num_diff_st)


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



x_k_mat<- matrix(numeric(920*920),920,920)
w_k_mat<- matrix(numeric(920*920),920,920)
x_w_k_mat<- matrix(numeric(920*920),920,920)

for(i in 1:920)
{
	for(j in 1:920)
	{
	x_k_mat[i,j]<- sum((all_x[i,]-all_x[j,])^2)
	w_k_mat[i,j]<- sum((all_w[i,]-all_w[j,])^2)
	x_w_k_mat[i,j]<- sum((all_x_w[i,]-all_x_w[j,])^2)

	}
}

length1<- median(x_k_mat)
length2<- median(w_k_mat)

sigmaf<- c(sigmaf1,sigmaf2)

compute_kernel<- function(a)
{
x1<- a[1:10]
w1<- a[11:51]
x2<- a[52:61]
w2<- a[62:102]


l<-a[103]

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

simple_PHV_regret_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))
simple_e_accuracy_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))
simple_e_accuracy2_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))

inference_PHV_regret_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))
inference_e_accuracy_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))
inference_e_accuracy2_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))

X_num_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))
AF_computation_time_mat<- matrix(numeric(total_seed_num*(iteration_num+50)),total_seed_num,(iteration_num+50))

###############################################################
###############################################################
###############################################################


###############################################################
##################     Optimization     #######################
###############################################################

Kvec_mat<-matrix(numeric(429*(iteration_num+50)),429,(iteration_num+50))%x%numeric(obj_dim)
post_mean_mat<- matrix(numeric(429*obj_dim),429,obj_dim)
post_var_mat<- matrix(numeric(429*obj_dim),429,obj_dim)
ucb_mat<- matrix(numeric(429*obj_dim),429,obj_dim)
lcb_mat<- matrix(numeric(429*obj_dim),429,obj_dim)
post_cov_mat<-diag(429)%x%numeric(obj_dim)

svd_mat<-diag(429)%x%numeric(obj_dim)

#original_all_x_w<-all_x_w

#all_x_w<-all_x_w*10


all_x_w<-all_x_w2[info[,2]==1,]
all_x_w<-scale(all_x_w)

#new_all_w<- matrix(numeric(429*41),429,41)
#for(i in 1:429)
#{
#new_all_w[i,]<-colMeans(rbind(all_w[info[,1]==i,]))
#
#
#}


#all_x<- scale(all_x)
#all_w<- scale(new_all_w)


x_w_k_mat<-diag(429)
#w_k_mat<-diag(429)

for(i in 1:429)
{
	for(j in 1:429)
	{
	x_w_k_mat[i,j]<- sum((all_x_w[i,]-all_x_w[j,])^2)
	#w_k_mat[i,j]<- sum((all_w[i,]-all_w[j,])^2)

	}
}

Len1<-c()
#Len2<-c()
for(i in 1:(total_input_num-1))
{
	Len1<-c(Len1,x_w_k_mat[i,(i+1):total_input_num])
	#Len2<-c(Len2,w_k_mat[i,(i+1):total_input_num])

}

length1<- median(Len1)
#length2<- median(Len2)



#all_x<- cbind(all_x,all_w)

all_x<- all_x_w


length_scale<-c(length1,length1)/2

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
new_X<- all_x[new_X_num,]
new_Y<- true_risk_mat[new_X_num,]#+rnorm(obj_dim,0,sqrt(sigman))
	
ww<-0
w4<-0

	while(w4<iteration_num)
	{
all_comp_start_time<-Sys.time()

	ww<-ww+1
	X<-rbind(X,new_X)	
	X_num<-c(X_num,new_X_num)
	Y<-rbind(Y,new_Y)

		for(www in 1:obj_dim)
		{
		Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),ww]<-t(matrix(apply(cbind(all_x,t(new_X)%x%(numeric(total_input_num)+1),www),1,compute_kernel),1,total_input_num))

		partial_Kvec_mat<-cbind(Kvec_mat[((www-1)*total_input_num+1):(www*total_input_num),1:ww])
		Cinv<-solve(partial_Kvec_mat[X_num,]+sigman*diag(ww))

		post_mean_mat[,www]<-as.vector(partial_Kvec_mat%*%Cinv%*%Y[,www])
		post_var_mat[,www]<- abs(sigmaf[www]-rowSums((partial_Kvec_mat%*%Cinv)*partial_Kvec_mat))

		}

	ucb_mat<- post_mean_mat+beta*sqrt(post_var_mat)
	lcb_mat<- post_mean_mat-beta*sqrt(post_var_mat)

	predicted_PF_PFn<-compute_PF(rbind(true_risk_mat[X_num,]))
	predicted_PF<-rbind(predicted_PF_PFn[,1:obj_dim])
	predicted_PFn<-X_num[predicted_PF_PFn[,(obj_dim+1)]]

predicted_PF_cell<-create_cells(predicted_PF,ref_point,apply(predicted_PF,2,max),c(1:length(predicted_PF[,1])))

	predicted_PF_PHV<-  sum(apply(rbind(predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-predicted_PF_cell[,(1):(obj_dim)]),1,prod))
	
	
	simple_PHV_regret<- c(simple_PHV_regret,max(true_PHV-predicted_PF_PHV,0))
	simple_PHV_regret_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-	min(simple_PHV_regret)

	simple_e_accuracy<- c(simple_e_accuracy,compute_e_accuracy(true_PF,predicted_PF))
	simple_e_accuracy_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-	min(simple_e_accuracy)

	simple_e_accuracy2<- c(simple_e_accuracy2,compute_e_accuracy2(true_PF,predicted_PF))
	simple_e_accuracy2_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-	min(simple_e_accuracy2)

	





	lcb_PF_num<-compute_PF(lcb_mat)[,(obj_dim+1)]
	reported_Y<-rbind(true_risk_mat[lcb_PF_num,])
	inference_predicted_PF_PFn<-compute_PF(rbind(true_risk_mat[lcb_PF_num,]))
	inference_predicted_PF<-rbind(inference_predicted_PF_PFn[,1:obj_dim])
	inference_predicted_PFn<-lcb_PF_num[inference_predicted_PF_PFn[,(obj_dim+1)]]

	inference_predicted_PF_cell<-create_cells(inference_predicted_PF,ref_point,apply(inference_predicted_PF,2,max),c(1:length(inference_predicted_PF[,1])))
	inference_predicted_PF_PHV<-  sum(apply(rbind(inference_predicted_PF_cell[,(obj_dim+1):(2*obj_dim)]-inference_predicted_PF_cell[,(1):(obj_dim)]),1,prod))


	inference_PHV_regret<- c(inference_PHV_regret,max(true_PHV-inference_predicted_PF_PHV,0))
	inference_PHV_regret_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-	inference_PHV_regret[ww]

	inference_e_accuracy<- c(inference_e_accuracy,compute_e_accuracy(true_PF,inference_predicted_PF))
	inference_e_accuracy_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-	inference_e_accuracy[ww]

	inference_e_accuracy2<- c(inference_e_accuracy2,compute_e_accuracy2(true_PF,reported_Y))
	inference_e_accuracy2_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-	inference_e_accuracy2[ww]

	X_num_mat[w,(w4+1):(w4+sum(info[,1]==new_X_num))]<-new_X_num
	w5<-w4
	w4<-(w4+sum(info[,1]==new_X_num))

	##########################################
		if(AF_num==1)
		{
		starttime<-Sys.time()

		new_X_num<- sample(c(1:total_input_num),1)

		endtime<-Sys.time()

		AF_computation_time_mat[w,(w5+1):(w5+sum(info[,1]==new_X_num))]<- difftime(endtime, starttime, units = "sec") 
		}

		if(AF_num==2)
		{
		starttime<-Sys.time()
		
		AFX<- apply(post_var_mat,1,sum)
		new_X_num<- order(AFX,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,(w5+1):(w5+sum(info[,1]==new_X_num))]<- difftime(endtime, starttime, units = "sec") 
		}

		if(AF_num==3)
		{
		starttime<-Sys.time()
		
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		current_PF<- rbind(rbind(compute_PF(Y))[,1:obj_dim])
		current_PF_cell<-create_cells(current_PF,apply(rbind(current_PF,ref_point),2,min),apply(current_PF,2,max),c(1:length(current_PF[,1])))
		current_PHV<- sum(apply(rbind(current_PF_cell[,(obj_dim+1):(2*obj_dim)]-current_PF_cell[,1:(obj_dim)]),1,prod))
		EMMDI_sample<- mvrnorm(EMMDI_sample_num,numeric(obj_dim),diag(obj_dim))

		AF_x<- function(a)
		{
		mu<- post_mean_mat[a,]
		sd<- sqrt(post_var_mat[a,])
		
		y1<- mu[1]+sd[1]*EMMDI_sample[,1]
		y2<- mu[2]+sd[2]*EMMDI_sample[,2]
		newY<- cbind(y1,y2)
		 compute_EMMDI(newY,current_PF)
		}


		AFX<- apply(cbind(c(1:429)),1,AF_x)

		new_X_num<- order(AFX,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,(w5+1):(w5+sum(info[,1]==new_X_num))]<- difftime(endtime, starttime, units = "sec") 
		}


		if(AF_num==4)
		{
		starttime<-Sys.time()
		
		set.seed(((seed-1)*total_seed_num+w)*1000+ww)
		
		LCB_PF<- rbind(rbind(compute_PF(lcb_mat))[,1:obj_dim])
		AFX1<-compute_dist(ucb_mat,LCB_PF)
		AFX<-  rowSums((ucb_mat-lcb_mat)^2)*(AFX1>=0)



		new_X_num<- order(AFX,decreasing=T)[1]

		endtime<-Sys.time()

		AF_computation_time_mat[w,(w5+1):(w5+sum(info[,1]==new_X_num))]<- difftime(endtime, starttime, units = "sec") 
		}

	##########################################

	new_X<- all_x[new_X_num,]
	new_Y<- true_risk_mat[new_X_num,]#+rnorm(obj_dim,0,sqrt(sigman))
		
all_comp_end_time<-Sys.time()
all_computation_time_mat[w,ww]<- difftime(all_comp_end_time, all_comp_start_time, units = "sec") 


	}

}



write.table(simple_PHV_regret_mat[,1:iteration_num],"simple_PHV_regret_mat.txt")
write.table(simple_e_accuracy_mat[,1:iteration_num],"simple_e_accuracy_mat.txt")
write.table(simple_e_accuracy2_mat[,1:iteration_num],"simple_e_accuracy2_mat.txt")


write.table(inference_PHV_regret_mat[,1:iteration_num],"inference_PHV_regret_mat.txt")
write.table(inference_e_accuracy_mat[,1:iteration_num],"inference_e_accuracy_mat.txt")
write.table(inference_e_accuracy2_mat[,1:iteration_num],"inference_e_accuracy2_mat.txt")


write.table(AF_computation_time_mat[,1:iteration_num],"AF_computation_time_mat.txt")
#write.table(X_num_mat,"X_num_mat.txt")

write.table(all_computation_time_mat,"all_computation_time_mat.txt")

