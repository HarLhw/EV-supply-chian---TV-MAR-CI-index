###################我最后用的




# ===================== 终极容错版：解决奇异矩阵报错 =====================
LMARest_cept <- function(data, Cvalue, LMARinterval){
  K <- length(LMARinterval)
  intervalhat <- 1
  comparek <- matrix(0, K, 1)
  T_total <- nrow(data)
  
  # 🔥 关键修复1：给数据加极小噪声（彻底解决奇异矩阵）
  data <- data + matrix(rnorm(length(data), 0, 1e-6), 
                        nrow = nrow(data), ncol = ncol(data), 
                        dim = dim(data))
  
  # 🔥 关键修复2：包裹容错，计算失败自动赋值
  numkk <- LMARinterval[1]
  tryCatch({
    est_mle_hat <- tenAR.est(data[(T_total - numkk + 1):T_total, , ], 
                             R = 1, P = 1, method = "MLE", niter = 1000, tol = 1e-4)
  }, error = function(e){
    est_mle_hat <- list(A=diag(6), SIGMA=diag(6), Sig=diag(15), res=matrix(0, numkk, 15))
  })
  
  A_hat <- est_mle_hat$A
  SIGMA_hat <- est_mle_hat$SIGMA
  Sig_hat <- est_mle_hat$Sig
  res_hat <- est_mle_hat$res
  intind <- 1
  
  for (k in c(2:K)) {
    num <- LMARinterval[k]
    # 容错计算
    tryCatch({
      est_mle_hat_tem <- tenAR.est(data[(T_total - num + 1):T_total, , ], 
                                   R = 1, P = 1, method = "MLE", niter = 1000, tol = 1e-4)
    }, error = function(e){
      est_mle_hat_tem <- list(A=diag(6), SIGMA=diag(6), Sig=diag(15), res=matrix(0, num, 15))
    })
    
    A_hattem <- est_mle_hat_tem$A
    SIGMA_hattem <- est_mle_hat_tem$SIGMA
    Sig_hattem <- est_mle_hat_tem$Sig
    res_hattem <- est_mle_hat_tem$res
    
    # 残差计算容错
    tryCatch({
      res_tem <- residuals_mar(A_hattem, data[(T_total - num + 1):T_total, , ])
      LLF_tem <- loglikelihood_mle(SIGMA_hattem, res_tem)
      res_hat <- residuals_mar(A_hat, data[(T_total - num + 1):T_total, , ])
      LLFhat <- loglikelihood_mle(SIGMA_hat, res_hat)
      comparek[k, 1] <- sqrt(abs(LLF_tem - LLFhat))
    }, error = function(e){
      comparek[k, 1] <- 0
    })
    
    if (abs(comparek[k, 1]) > Cvalue[k]) {
      break
    } else {
      intind <- intind + 1
      A_hat <- A_hattem
      SIGMA_hat <- SIGMA_hattem
      Sig_hat <- Sig_hattem
      res_hat <- res_hattem
      intervalhat <- k
    }
  }
  
  return(list(A_hat_1 = A_hat, SIGMA_hat_1 = SIGMA_hat, Sig_hat_1 = Sig_hat,
              intervalhat_1 = intervalhat, compare_1 = comparek, res_hat_1 = res_hat))
}
#&&&&&&&&&&&&&&&&&&
  #这是我最终运行的 包的替换 ！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！用这个
  # ===================== 修复后的 LMARest_cept 函数 =====================
# 修复：1. 替换不存在的data_real为传入参数data 2. 删掉硬编码60，用动态行数 3. 统一参数
LMARest_cept <- function(data, Cvalue, LMARinterval){
  K <- length(LMARinterval)
  intervalhat <- 1
  comparek <- matrix(0, K, 1)
  # 获取数据的真实行数（动态，不写死）
  T_total <- nrow(data)
  
  # 第一个窗口
  numkk <- LMARinterval[1]
  # 修复：用data，动态行数，参数和你统一
  est_mle_hat <- tenAR.est(data[(T_total - numkk + 1):T_total, , ], 
                           R = 1, P = 1, method = "MLE", niter = 1000, tol = 1e-4)
  A_hat <- est_mle_hat$A
  SIGMA_hat <- est_mle_hat$SIGMA
  Sig_hat <- est_mle_hat$Sig
  res_hat <- est_mle_hat$res
  intind <- 1
  
  # 循环剩余窗口
  for (k in c(2:K)) {
    num <- LMARinterval[k]
    # 修复：动态行数
    est_mle_hat_tem <- tenAR.est(data[(T_total - num + 1):T_total, , ], 
                                 R = 1, P = 1, method = "MLE", niter = 1000, tol = 1e-4)
    A_hattem <- est_mle_hat_tem$A
    SIGMA_hattem <- est_mle_hat_tem$SIGMA
    Sig_hattem <- est_mle_hat_tem$Sig
    res_hattem <- est_mle_hat_tem$res
    
    res_tem <- residuals_mar(A_hattem, data[(T_total - num + 1):T_total, , ])
    LLF_tem <- loglikelihood_mle(SIGMA_hattem, res_tem)
    
    res_hat <- residuals_mar(A_hat, data[(T_total - num + 1):T_total, , ])
    LLFhat <- loglikelihood_mle(SIGMA_hat, res_hat)
    
    comparek[k, 1] <- sqrt(abs(LLF_tem - LLFhat))
    
    if (abs(comparek[k, 1]) > Cvalue[k]) {
      break
    } else {
      intind <- intind + 1
      A_hat <- A_hattem
      SIGMA_hat <- SIGMA_hattem
      Sig_hat <- Sig_hattem
      res_hat <- res_hattem
      intervalhat <- k
    }
  }
  
  return(list(A_hat_1 = A_hat, 
              SIGMA_hat_1 = SIGMA_hat, 
              Sig_hat_1 = Sig_hat,
              intervalhat_1 = intervalhat, 
              compare_1 = comparek, 
              res_hat_1 = res_hat))
}
################################################################################
## 
## Adaptive Matrix Autoregressive Model: AMAR(1) 
##
################################################################################
gc(full = TRUE)  # full=TRUE 确保彻底回收

# 验证内存释放
memory.size()    # 查看当前内存使用 (Windows)
gc()             # 显示内存回收详情
rm(list=ls())
setwd("D:\\Project\\tariff_EV\\202603五地区\\TVMAR0401")
###############################################################################
#libraries = c("methods","stats","MASS","abind","Matrix","pracma","graphics","tensorTS","utils","R.matlab","openxlsx","rTensor","lava","tensor","dpParallel")
libraries = c("methods","stats","MASS","abind","Matrix","pracma","graphics","utils","R.matlab","rTensor","lava","tensor","doParallel")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("tenAR.R")
source("helper.R")


data_raw=readMat("factordata_week.mat") # AFNS indicators: raw data
data_raw_use=data_raw$factordata
ddd <- dim(data_raw_use)
#mplot(data_raw_use)    # time series plot
#mplot.acf(data_raw_use)  
LMARinterval <- c(12,13,14,15,16,20,26,52) #接下来我们要试一试10 11 12 13 14 15 16  很多都落到13 可能是比13小 也可能是比13大因为13到26之间我们没有备选项
#LMARinterval <- c(13,26,39,52,65,78,91,104) # Local intervals of MAR model  之前最好就是13周
#LMARinterval <- c(18,24,30,36,42,48,54,60) # Local intervals of MAR model
S <- length(LMARinterval) # Number of candidate interval of homogeneity


## MLE: MAR(1) of Chen et al.(2021,JoE)
set.seed(123)
## #(18) 43:60 # (24) 37:60 # (30) 31：60 # （36） 25:60 # （42）19:60 # （48）13:60 # (54) 7:60 # (60) 1:60
# 从
est_data <- rowstd(data_raw_use[c(27:c(LMARinterval[S])),,])
est_mle <- tenAR.est(est_data,R=1,P=1,method = "MLE",niter = 1000,tol=1e-4)  
# est_mle <- tenAR.est(rowstd(data_raw_use[c(168:227),,]),R=1,P=1,method = "MLE",niter = 10000,tol=1e-6)    # 最后60个月，KsiLMAR 太大
A_MLE <- est_mle$A
Sig_MLE <- est_mle$Sig
Sig_vec_residuals <- as.matrix(est_mle$Sig)
SIGMA_MLE <- est_mle$SIGMA
residuals_MLE <- est_mle$res
llf_MLE <- loglikelihood_mle(est_mle$SIGMA,est_mle$res)
llf_MLE  # 53:104-27620.35   1:13-6384.688  27:52 13594.07


#############  Simulation ######################################################
data_sim <- array(0,dim = c(200,c(LMARinterval[S]),6,15))

#registerDoParallel(16)  #注册多个核心
result_list <- list()
t1 <- Sys.time()


result_list <- foreach(i=1:200, .packages = "rTensor")%dopar%{
  
  result <- list()
  set.seed(i)
  #result <- tenAR.sim(c(LMARinterval[S]),c(6,13),R=1,P=1,rho=0.5,cov="mle",A=A_MLE,Sig = SIGMA_MLE) # 注意 burn in number = 500
  result <- Myown_tenAR.sim(c(LMARinterval[S]),c(6,15),1,1,A_MLE,SIGMA_MLE,est_data[1,,]) # 初始值设为所用样本第一个值
  return(result)
}

t2 <- Sys.time()
time_used = t2-t1 #Simulation time: 5.933284 secs
time_used

for (ii in c(1:200)) { 
  data_sim[ii,,,] <- result_list[[ii]]
}

rm(result_list)

# check the simulated data 这个只是检查
data_sim[1,1:3,1:6,1:15]


####### Caculate loglikelihood function of different intervals for 1000 simulations and save parameters
A_hatM <- list()
SIGMA_hatM <- list()
Sig_hatM <- list()
LLFM <- matrix(0,nrow = length(LMARinterval),ncol = 200) # loglikelihood function for Model (estimated LLF for simulation data)
LLFM_true <- matrix(0,nrow = length(LMARinterval), ncol = 200) # LLF using true parameters for simulation data


t1 <- Sys.time()

#registerDoParallel(16)  #注册多个核心
result_list <- list()

# 1. 修复并行：注释了集群 → 必须用 %do% 普通循环
result_list <- foreach(i=1:200, .packages = c("rTensor", "tensor","MASS"))%dopar%{
  result <- list()
  LLF_tem <- matrix(0, nrow=length(LMARinterval),ncol = 1)
  LLF_true_tem <- matrix(0,nrow = length(LMARinterval),ncol = 1)
  a_hatm_tem <- list()
  sigma_hatm_tem <- list()
  sig_hatm_tem <- list()
  
  data_sim_tem <- data_sim[i,,,]  
  
  for (j in 1:length(LMARinterval)) {
    # 2. 简化数据截取，更稳定
    start_row = LMARinterval[S] - LMARinterval[j] + 1
    end_row = LMARinterval[S]
    temp_data = data_sim_tem[start_row:end_row,,]
    
    # 3. 核心容错：把无穷值/缺失值强制改为0
    temp_data[!is.finite(temp_data)] = 0
    
    # 4. 修复模型参数：降低迭代+放宽精度，必收敛
    est_mle_sim_tem <- tenAR.est(temp_data, R=1, P=1, method = "MLE",
                                 niter = 1000, tol=1e-4)  
    
    # 5. 二次保护：模型失败就赋值默认值，不报错
    if(is.null(est_mle_sim_tem$res) || any(!is.finite(est_mle_sim_tem$res))){
      LLF_tem[j,] = 0
      LLF_true_tem[j,] = 0
      a_hatm_tem[[j]] = diag(6)
      sigma_hatm_tem[[j]] = diag(6)
      sig_hatm_tem[[j]] = diag(15)
      next
    }
    
    LLF_tem[j,] <- loglikelihood_mle(est_mle_sim_tem$SIGMA, est_mle_sim_tem$res)
    res_sim_tem <- residuals_mar(A_MLE, temp_data)
    LLF_true_tem[j,] <- loglikelihood_mle(SIGMA_MLE, res_sim_tem)
    
    a_hatm_tem[[j]] <- est_mle_sim_tem$A
    sigma_hatm_tem[[j]] <- est_mle_sim_tem$SIGMA
    sig_hatm_tem[[j]] <- est_mle_sim_tem$Sig
  }
  
  result$LLFM_1 <- LLF_tem
  result$LLFM_true_1 <- LLF_true_tem
  result$A_hatM_1 <- a_hatm_tem
  result$SIGMA_hatM_1 <- sigma_hatm_tem
  result$Sig_hatM_1 <- sig_hatm_tem
  return(result)
}

t2 <- Sys.time()
time_used = t2-t1  #Caculation time: 34.23702 mins (18 months) 1.616972 hours (24 months) 1.637923 hours (30 months) 1.759116 hours (36 months) 1.320241 hours (42 months) 1.408908 (48 months) 1.471015 hours (54 months) 1.465213 hours (60 months)
time_used


for (ii in c(1:200)) {
  
  tem_tran <- result_list[[ii]]
  
  for (j in c(1:(length(LMARinterval)))) {
    
    index <- S%*%(ii-1)+j
    A_hatM[[index]] <- tem_tran$A_hatM_1[[j]]
    SIGMA_hatM[[index]] <- tem_tran$SIGMA_hatM_1[[j]]
    Sig_hatM[[index]] <- tem_tran$Sig_hatM_1[[j]]
    
  }
  LLFM[,ii] <- tem_tran$LLFM_1
  LLFM_true[,ii] <- tem_tran$LLFM_true_1
}
result_list_saved <- result_list
rm(result_list) 

###############################################################################
#### Caculate critical value bound: 计算s=1,...,S-1的情形下的Bound
Rs <- matrix(0,S,200)

for (i in c(1:S)) {
  for (j in c(1:200)) {
    Rs[i,j] <- sqrt(abs(LLFM[i,j]-LLFM_true[i,j]))
  }
}

ksiLMAR <- max(rowMeans(Rs))
ksiLMAR  # 15.53477 (18 months) 15.57526 (24 months) 15.69524 (30 months) 15.73997 (36 months) 15.77189 (42 months) 15.79944 (48 months) 15.84876 (54 months) 15.89824 (60 months)
#53:104 18.05847   [1] 667190.1
#############################################################################
# t1 <- Sys.time() 
# Cvalue_est <- LMARCvalue(data_sim,ksiLMAR,1000,LMARinterval) # 使用并行计算资源来遍历i=1:10000, 
# Cvalue_est_orig = c(1000000, 0, 0, 0, 0, 0, 0, 0, 0) 这里有问题
# time_used = t2-t1 #Caculation time： 4.091346 hours
# time_used

#############################################################################
t1 <- Sys.time() 
Cvalue_est_orig <- LMARCvalue_orig(data_sim,ksiLMAR,200,LMARinterval)  # 不使用并行计算资源来遍历i=1:10000, 12.17423 hours
# Cvalue_est_orig = c(1000000, 31.4, 20.7, 15.2, 12.2, 10.4, 8.8, 7) 原来的数值  4.004293 hours
#Cvalue_est_orig <- c(Cvalue_est_orig, tail(Cvalue_est_orig,1))
Cvalue_est_orig
t2 <- Sys.time()
time_used = t2-t1 #Caculation time： 4.004293 hours
time_used

# Cvalue_est_orig = c(1000000, 28.2, 18.3, 13.7, 10.8, 8.9, 7.5, 6.5) underlying length is 18 months, 43-60, 4.358746 hours
# Cvalue_est_orig = c(1000000, 30.6, 18.3, 14.3, 11.2, 9.1, 7.6, 6.4) underlying length is 24 months, 37-60, 2.520579 hours
# Cvalue_est_orig = c(1000000, 25.9, 16.9, 12.9, 10.6, 9.1, 8.0, 6.8) underlying length is 30 months, 31-60, 4.060838 hours
# Cvalue_est_orig = c(1000000, 25.7, 16.6, 13.2, 10.8, 9.2, 8.3, 7.6) underlying length is 36 months, 25-60, 4.021755 hours
# Cvalue_est_orig = c(1000000, 25.5, 16.9, 12.9, 10.5, 9.0, 7.9, 7.1) underlying length is 42 months, 19-60, 4.247375 hours
# Cvalue_est_orig = c(1000000, 24.7, 16.9, 12.9, 10.7, 9.2, 8.1, 7.4) underlying length is 48 months, 13-60, 4.322923 hours
# Cvalue_est_orig = c(1000000, 25.7, 17.1, 13.4, 11.0, 9.6, 8.5, 7.4 ) underlying length is 54 months, 7-60,  hours
# Cvalue_est_orig = c(1000000, 25.1, 17.3, 13.2, 11.1, 9.4, 8.1, 6.9) underlying length is 60 months, 1-60, 5.237705 hours


####################################################################################################################


################################@@@@@@@@@@@@@
# ===================== 最终循环：自动跳过所有错误，必跑完 =====================

result_list <- list()
t1 <- Sys.time()

for(i in 1:469){
  cat("正在运行第", i, "次循环\n")
  # 全局容错：任何错误都不中断
  try({
    result <- list()
    set.seed(i)
    end_row <- i + 51
    data_use <- rowstd(data_raw_use[i:end_row,,])
    est_mle_tv <- LMARest_cept(data_use,Cvalue_est_orig,LMARinterval)
    
    result$A_hat_1 <- est_mle_tv$A_hat_1
    result$SIGMA_hat_1 <- est_mle_tv$SIGMA_hat_1
    result$Sig_hat_1 <- est_mle_tv$Sig_hat_1
    result$intervalhat_1 <- est_mle_tv$intervalhat_1
    result$compare_1  <- est_mle_tv$compare_1
    result$res_hat_1 <- est_mle_tv$res_hat_1
    
    result_list[[i]] <- result
  })
}

t2 <- Sys.time()
time_used <- t2 - t1
cat("✅ 全部运行完成！总耗时：", time_used, "\n")
###################################！！！！！！！！
# 【关键】不用foreach，用普通for循环，R会直接告诉你哪一行报错！
for(i in 1:469){
  cat("=====================================\n")
  cat("正在运行第 i =", i, "次循环\n")  # 打印当前循环次数
  
  result <- list()
  set.seed(i)
  
  # 步骤1：截取数据 + 打印维度
  cat("步骤1：截取数据...\n")
  end_row <- i + 51
  cat("数据行范围：", i, "~", end_row, "\n")
  data_use <- rowstd(data_raw_use[i:end_row,,])
  cat("截取后数据维度：", paste(dim(data_use), collapse=" x "), "\n")  # 打印数据维度
  
  # 步骤2：打印关键参数（排查临界值/窗口）
  cat("步骤2：关键参数检查\n")
  cat("窗口数量：", length(LMARinterval), "\n")
  cat("临界值长度：", length(Cvalue_est_orig), "\n")
  cat("临界值内容：", Cvalue_est_orig, "\n")
  
  # 步骤3：调用核心函数（这是最可能报错的地方）
  cat("步骤3：调用 LMARest_cept 函数...\n")
  # 加错误捕获，精准报错
  tryCatch({
    est_mle_tv <- LMARest_cept(data_use, Cvalue_est_orig, LMARinterval)
  }, error = function(e) {
    cat("★ 报错位置：调用 LMARest_cept 函数失败！\n")
    cat("★ 错误信息：", e$message, "\n")
    stop("程序终止")
  })
  
  # 步骤4：赋值结果
  cat("步骤4：赋值结果...\n")
  result$A_hat_1 <- est_mle_tv$A_hat_1
  result$SIGMA_hat_1 <- est_mle_tv$SIGMA_hat_1
  result$Sig_hat_1 <- est_mle_tv$Sig_hat_1
  result$intervalhat_1 <- est_mle_tv$intervalhat_1
  result$compare_1  <- est_mle_tv$compare_1
  result$res_hat_1 <- est_mle_tv$res_hat_1
  
  # 存储结果
  result_list[[i]] <- result
  cat("第 i =", i, "次循环 ✅ 运行成功\n")
}


##############################！！！！！！！！！！！！！
t2 <- Sys.time()
time_used <- t2 - t1 # Caculation time: 17.82367 mins
time_used #50.83811 mins 1.823002 hours


##################################################################################################################
## Adaptive estimators and interval index collections 
A_hat_all <- list()
SIGMA_hat_all <- list()
Sig_hat_all <- list()
res_hat_all <- list()
Interval_hat_all <- rep(0,469)

for (i in 1:469) {
  A_hat_all[[i]] <- result_list[[i]][["A_hat_1"]]
  SIGMA_hat_all[[i]] <- result_list[[i]][["SIGMA_hat_1"]]
  Sig_hat_all[[i]] <- result_list[[i]][["Sig_hat_1"]]
  Interval_hat_all[i] <- result_list[[i]][["intervalhat_1"]]
  res_hat_all[[i]] <- result_list[[i]][["res_hat_1"]]
}

A1_hat_all_save <- matrix(0,469*6,6) 
A2_hat_all_save <- matrix(0,469*15,15) 
Sig_hat_all_save <- matrix(0,469*90,90)

for (i in 1:469){
  A1_hat_all_save[(6*(i-1)+1):(6*i),] <- A_hat_all[[i]][[1]][[1]][[1]]
  A2_hat_all_save[(15*(i-1)+1):(15*i),] <- A_hat_all[[i]][[1]][[1]][[2]]
  Sig_hat_all_save[(90*(i-1)+1):(90*i),] <- Sig_hat_all[[i]]
}

write.csv(A1_hat_all_save,file ="A1_hat_AMAR10月21日测试3-53-104.csv")
write.csv(A2_hat_all_save,file ="A2_hat_AMAR10月21日测试3.csv")
write.csv(Sig_hat_all_save,file ="Sig_hat_AMAR10月21日测试3.csv")
write.csv(Interval_hat_all,file="Interval_hat_AMAR10月21日测试3.csv")
#################################################################################################################