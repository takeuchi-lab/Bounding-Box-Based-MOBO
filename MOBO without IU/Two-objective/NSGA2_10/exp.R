library(MASS)
library(nsga2R)
seed<-1
total_seed_num<-100
iteration_num<-30
PHV_sample_num<-10000
AF_num<-5			#1:Random, 2:US, 3:EHVI, 4:e-PAL, 5:Proposed, 6:modified_proposed, 7:ParEGO, 8:SMSego, 9:MESMO, 10:PFES


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

###############################################################
##################     True function    #######################
###############################################################

true_function<-function(a)
{
x<-a[1]
y<-a[2]
	
	f1<- (-1*( (x+2*y-7)^2+(2*x+y-5)^2 )-(-157.35))/sqrt(28896.11)
	


	f2<- (-1*( 0.26*(x^2+y^2)-0.48*x*y )-(-4.3342))/sqrt(23.52052)
	c(f1,f2)
}


negative_true_function<-function(a)
{
x<-a[1]
y<-a[2]
	
	f1<- (-1)*(-1*( (x+2*y-7)^2+(2*x+y-5)^2 )-(-157.35))/sqrt(28896.11)
	


	f2<- (-1)*(-1*( 0.26*(x^2+y^2)-0.48*x*y )-(-4.3342))/sqrt(23.52052)
	c(f1,f2)
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

true_function_mat<- t(apply(all_input,1,true_function))

ref_point<-apply(true_function_mat,2,min)
max_point<-apply(true_function_mat,2,max)
vol_refp_maxp<- prod(max_point-ref_point)

###############################################################
###############################################################
###############################################################


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





#####################################################

true_PF_PFn<-compute_PF(true_function_mat)
true_PF<-true_PF_PFn[,1:obj_dim]
true_PFn<-true_PF_PFn[,(obj_dim+1)]

#set.seed(1)
true_PF_cell<- create_cells(true_PF,ref_point,max_point,c(1:length(true_PF[,1])))
true_PHV<- sum(apply(rbind(rbind(true_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(true_PF_cell)[,1:(obj_dim)]),1,prod))


###############################################################
##################      Result mat      #######################
###############################################################

simple_PHV_regret_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
simple_e_accuracy_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
simple_e_accuracy2_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)

inference_PHV_regret_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
inference_e_accuracy_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)
inference_e_accuracy2_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)


AF_computation_time_mat<- matrix(numeric(total_seed_num*iteration_num),total_seed_num,iteration_num)





for(w in 1:total_seed_num)
{


	for(ww in 1:iteration_num)
	{
	set.seed(w)

	starttime<-Sys.time()

	

	results <- nsga2R(fn=negative_true_function, varNo=2, objDim=2, lowerBounds=c(lower,lower), upperBounds=c(upper,upper),
	popSize=10, tourSize=2, generations=ww, cprob=0.9, XoverDistIdx=20, mprob=0.1,MuDistIdx=3)
	reported_Y<- t(apply((results$parameters),1,true_function))
	NSGA2_PF<- rbind(compute_PF(reported_Y)[,1:obj_dim])

	NSGA2_PF_cell<- create_cells(NSGA2_PF,ref_point,max_point,c(1:length(NSGA2_PF[,1])))
	NSGA2_PF_PHV<- sum(apply(rbind(rbind(NSGA2_PF_cell)[,(obj_dim+1):(2*obj_dim)]-rbind(NSGA2_PF_cell)[,1:(obj_dim)]),1,prod))







	

	inference_PHV_regret_mat[w,ww]<-	max(true_PHV-NSGA2_PF_PHV,0)
	simple_PHV_regret_mat[w,ww]<-	min(inference_PHV_regret_mat[w,1:ww])


	inference_e_accuracy_mat[w,ww]<-	compute_e_accuracy(true_PF,NSGA2_PF)


	inference_e_accuracy2_mat[w,ww]<-	compute_e_accuracy2(true_PF,reported_Y)

	
	endtime<-Sys.time()

	AF_computation_time_mat[w,ww]<- difftime(endtime, starttime, units = "sec") 

	}


}


write.table(simple_PHV_regret_mat,"simple_PHV_regret_mat.txt")
write.table(simple_e_accuracy_mat,"simple_e_accuracy_mat.txt")
write.table(simple_e_accuracy2_mat,"simple_e_accuracy2_mat.txt")


write.table(inference_PHV_regret_mat,"inference_PHV_regret_mat.txt")
write.table(inference_e_accuracy_mat,"inference_e_accuracy_mat.txt")
write.table(inference_e_accuracy2_mat,"inference_e_accuracy2_mat.txt")


write.table(AF_computation_time_mat,"AF_computation_time_mat.txt")





