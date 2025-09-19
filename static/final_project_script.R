#load necessary packages
library(ggplot2)
library(reshape2)
library(dplyr)
library(here)
library(socialmixr)
library(odin)

data<-read.csv(here("data","daily_GE_data.csv"))
data$t<-as.factor(data$Days)
# Barplot
ggplot(data, aes(x=t, y=Cases)) + 
  geom_bar(stat = "identity", width=0.2)+
  ylab("reported GE cases") +
  xlab("day")

## 2. NoV bulaşması için salgın verilerine bir model kalibre edin

## NoV modeli oluştur 
seiar_generator <- odin::odin({
  
  dt <- user(1)
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  
  ## Core equations for overall transitions between compartments:
  update(V_tot) <- V_tot + sum(n_SV) + sum(n_RV) - sum(n_VS) - sum(n_VVA) 
  update(S_tot) <- S_tot + sum(n_VS) + sum(n_RS) - sum(n_SE) - sum(n_SV)  
  update(E_tot) <- E_tot + sum(n_SE)  - sum(n_EI) 
  update(I_tot) <- I_tot + sum(n_EI) - sum(n_IA) - sum(n_deathI)
  update(A_tot) <- A_tot + sum(n_IA) - sum(n_AR) - sum(n_deathA)
  update(VA_tot) <- VA_tot +  sum(n_VVA) - sum(n_VAR) - sum(n_deathVA)
  update(R_tot) <- R_tot + sum(n_AR) + sum(n_VAR) - sum(n_RS) - sum(n_RV) 
  
  ## Equations for transitions between compartments by age group
  update(V[]) <- V[i] + n_SV[i] + n_RV[i] - n_VS[i] - n_VVA[i] 
  update(S[]) <- S[i] + n_VS[i] + n_RS[i] - n_SE[i] - n_SV[i]    
  update(E[]) <- E[i] + n_SE[i] - n_EI[i] 
  update(I[]) <- I[i] + n_EI[i] - n_IA[i] - n_deathI[i]
  update(A[]) <- A[i] + n_IA[i] - n_AR[i] - n_deathA[i]
  update(VA[]) <- VA[i] + n_VVA[i] - n_VAR[i] - n_deathVA[i]
  update(R[]) <- R[i] + n_AR[i] + n_VAR[i] - n_RS[i] - n_RV[i]
  
  ## Model outputs
  update(cum_vaccines[]) <-cum_vaccines[i]+ n_SV[i] + n_RV[i]
  update(new_cases[]) <- n_EI[i] 
  update(cum_cases[]) <- cum_cases[i] + n_EI[i] 
  update(cum_deaths[]) <- cum_deaths[i] + n_deathI[i] + n_deathA[i] + n_deathVA[i]
  
  update(cum_vaccines_all) <- cum_vaccines_all+ sum(n_SV) + sum(n_RV)
  update(new_cases_all) <- sum(n_EI)
  update(cum_cases_all) <- cum_cases_all + sum(n_EI) 
  update(new_reported_all) <- sum(new_reports)
  update(cum_deaths_all) <- cum_deaths_all + sum(n_deathI) + sum(n_deathA) + sum(n_deathVA)
  
  new_reports[]<- n_EI[i] * 1/rep_ratio[i]
  dim(new_reports)<-N_age
  
  
  ## Individual probabilities of transition:
  p_VS[] <- 1 - exp(-delta * dt)  # V to S
  p_SE[] <- 1 - exp(-lambda[i] * dt) # S to E
  p_EI   <- 1 - exp(-epsilon * dt) # E to I
  p_IA   <- 1 - exp(-theta * dt) # I to A
  p_AR   <- 1 - exp(-sigma * dt) # A to R
  p_RS   <- 1 - exp(-tau * dt) # R to S
  p_vacc[] <- 1 - exp(-vac_imm*vac_cov[i] * dt)# vaccination
  p_noromu[]<- 1 - exp(- ((p_IA+p_AR) * cfr[i]/(1-cfr[i])) * dt)# vaccination
  
  ## Force of infection
  m[, ] <- user() # age-structured contact matrix
  s_ij[, ] <- m[i, j] * (I[j] + (A[j] * rho )+ (A[j] * rho * vac_eff) )
  lambda[] <- beta * sum(s_ij[i, ])
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  
  # Flowing out of V 
  n_VS[] <- rbinom(V[i],  p_VS[i])
  n_VVA[] <- rbinom(V[i]-n_VS[i], p_SE[i] * (1-vac_eff))
  
  # Flowing out of S
  n_SE[] <- rbinom(S[i], p_SE[i])
  n_SV[] <- rbinom(S[i]-n_SE[i], if (step > t_vacc) p_vacc[i] else 0)
  
  # Flowing out of E
  n_EI[] <- rbinom(E[i], p_EI)
  
  # Flowing out of I
  n_IA[] <- rbinom(I[i], p_IA)
  n_deathI[]<-rbinom(I[i] - n_IA[i], p_noromu[i])
  
  # Flowing out of A
  n_AR[] <- rbinom(A[i], p_AR)
  n_deathA[]<-rbinom(A[i] - n_AR[i], p_noromu[i])
  
  # Flowing out of VA
  n_VAR[] <- rbinom(VA[i], p_AR)
  n_deathVA[]<-rbinom(VA[i] - n_VAR[i], p_noromu[i]* (1-vac_eff))
  
  # Flowing out of R
  n_RS[] <- rbinom(R[i], p_RS)
  n_RV[] <- rbinom(R[i] - n_RS[i], if (step > t_vacc) p_vacc[i] else 0)
  
  ## Initial states:
  initial(V_tot) <- sum(V_ini)
  initial(S_tot) <- sum(S_ini)
  initial(E_tot) <- sum(E_ini)
  initial(I_tot) <- sum(I_ini)
  initial(A_tot) <- sum(A_ini)
  initial(VA_tot) <- sum(VA_ini)
  initial(R_tot) <- sum(R_ini)
  
  initial(V[]) <- V_ini[i]
  initial(S[]) <- S_ini[i]
  initial(E[]) <- E_ini[i]
  initial(I[]) <- I_ini[i]
  initial(A[]) <- A_ini[i]
  initial(VA[]) <- VA_ini[i]
  initial(R[]) <- R_ini[i]
  initial(cum_vaccines[]) <- 0
  initial(new_cases[]) <- 0
  initial(cum_cases[]) <- 0
  initial(cum_deaths[]) <- 0
  initial(cum_vaccines_all) <- 0
  initial(new_cases_all) <- 0
  initial(cum_cases_all) <- 0
  initial(new_reported_all)<-0
  initial(cum_deaths_all) <- 0
  
  
  ## User defined parameters - default in parentheses:
  V_ini[] <- user()
  S_ini[] <- user()
  E_ini[] <- user()
  I_ini[] <- user()
  A_ini[] <- user()
  VA_ini[] <- user()
  R_ini[] <- user()
  
  # Model parameters (values in brackets are default values) 
  beta <- user(0.003)   # transm coefficient
  delta <- user( 1/(365))  # vaccine immunity dur (~5 yras)
  epsilon <- user(1)   # incubation
  theta <- user(0.5)   # duration symptoms
  sigma <- user(0.066) # duration asymp shedding
  tau   <- user(1/365)  # duration immunity
  rho   <- user(0.05) # rel infect asymptomatic 
  cfr[]  <- user() # Noro CFR by age 
  vac_eff<-user(0.9) # Vaccine efficacy for transmission
  vac_cov[]<-user() # vaccine coverage by age group
  vac_imm  <- user(1/5)# time to vaccine seroconversion (days)
  t_vacc   <- user(2) # days after case 0 to intro vaccine
  rep_ratio[]  <-user() # cases in community per case reported 
  # dimensions of arrays
  N_age <- user()
  dim(V_ini) <- N_age
  dim(S_ini) <- N_age
  dim(E_ini) <- N_age
  dim(I_ini) <- N_age
  dim(A_ini) <- N_age
  dim(VA_ini) <- N_age
  dim(R_ini) <- N_age
  dim(vac_cov)<-N_age
  dim(cfr)  <- N_age  
  dim(rep_ratio)  <- N_age  
  dim(V) <- N_age
  dim(S) <- N_age
  dim(E) <- N_age
  dim(I) <- N_age
  dim(A) <- N_age
  dim(VA) <- N_age
  dim(R) <- N_age
  dim(cum_vaccines) <- N_age
  dim(new_cases) <- N_age
  dim(cum_cases) <- N_age
  dim(cum_deaths) <- N_age
  dim(n_SV) <- N_age
  dim(n_RV) <- N_age
  dim(n_VS) <- N_age
  dim(n_VVA) <- N_age
  dim(n_VAR) <- N_age
  dim(n_SE) <- N_age
  dim(n_EI) <- N_age
  dim(n_IA) <- N_age
  dim(n_AR) <- N_age
  dim(n_RS) <- N_age
  dim(n_deathI)<-N_age
  dim(n_deathA)<-N_age
  dim(n_deathVA)<-N_age
  dim(p_VS) <- N_age
  dim(p_SE) <- N_age
  dim(p_vacc)<-N_age
  dim(p_noromu)<-N_age
  dim(m) <- c(N_age, N_age)
  dim(s_ij) <- c(N_age, N_age)
  dim(lambda) <- N_age
}, verbose = FALSE)

# Temas matrisini yükle 
cmat<-read.csv(here("data","contact_TUR.csv"))

contact_matrix<-as.matrix(cmat) # Matris nesnesine dönüştür
rownames(contact_matrix)<-c("0-4","5-14","14-64","65+") # yaş etiketleri ekle
colnames(contact_matrix)<-c("0-4","5-14","14-64","65+") # yaş etiketleri ekle

# temas matrisini çizmek için matrix_plot kullan
matrix_plot(contact_matrix)


# Rastgele sayılar için tohum
set.seed(1)

# İlgilenilen X bölgesinde toplam nüfusu tanımla
N<-68000

n_age<- 4 # yaş grubu sayısı

# Nüfus parametreleri
#  0-4  5-14  15-64  65+
pop_distr<- c(0.16, 0.17, 0.63,  0.04) # Nüfus yaş dağılımı
pop <- round(N * pop_distr)

# Temas matrisini işle: matrisin simetrik olduğundan emin olmamız gerekiyor 
# yani talep edilen temaslar sunulan temaslar ile eşit
cmat_sym<-((cmat+t(cmat))/2)

# Bulaş formülüne girdi olacak kişi başına temas oranını bul 
# yani her gruptaki nüfus büyüklüğüne göre düzeltildi
transmission <- as.matrix(cmat_sym )/
  rep(c(t(pop)), each = ncol(cmat_sym))

# Noro Vaka ölüm oranı (bkz. Lindsay et al https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-015-1168-5))

# Yaş grubuna göre NoV CFR
#  0-4  5-14  15-64  65+
mu<-c( 0.04, 0.01, 0.03, 0.63)/1000

# Bildirilen oranlar: hastane salgınında bildirilen vaka başına toplumda vaka sayısı
# (Varsayım: bu kesin değil ama örn. Birleşik Krallık'ta bildirilen her vaka için 
# toplumda ~280 bulunabileceği tahmin ediliyor)

#  0-4  5-14  15-64  65+
rep_ratio<-c(40,  65,   30,    15)

# SEIAR nesnesini çağır
seiar0 <- seiar_generator$new(
  V_ini = c(1:n_age)*0,
  S_ini = as.numeric(round(N*pop_distr - c(1,0,0,0))),
  E_ini = c(1:n_age)*0,
  I_ini = c(1,0,0,0),
  A_ini = c(1:n_age)*0,
  VA_ini = c(1:n_age)*0,
  R_ini = c(1:n_age)*0,
  N_age = n_age,
  cfr= mu,
  m=transmission,
  rep_ratio=rep_ratio,
  vac_cov= c(0,0,0,0),
  beta = 1 ## <================ Aşağıdaki grafikte en iyi uyan farklı beta değerlerini deneyin
)


# Çoklu çalıştırmalar (100)
t_end<- 365 * 2 # sim zamanı (2 yıl)

# Modeli çalıştır
seiar0_100 <- seiar0$run(0:t_end, replicate = 100)

# Değişken indeksi
idx<-rownames(seiar0_100[1,,])

# Bildirilen vakalar vs veri : Doğru betayı bulana kadar bu kodda yineleme yapın
t_id<-which(idx=="time")
id<- which(idx=="new_reported_all")
mean <- rowMeans(seiar0_100[, id,])
matplot(seiar0_100[, t_id,],seiar0_100[, id,], 
        xlab = "Günler", 
        ylab = "Bildirilen GE vaka sayısı",
        type = "l", lty = 1, col="grey",
        xlim=c(0,30),
        ylim=c(0,max(data$Cases)*1.2))
lines(seiar0_100[, 1,1],mean,col="purple")
points(data$Days+5, data$Cases, col = "red", pch = 19)

# Toplum insidansı
id<- which(idx=="new_cases_all")
mean <- rowMeans(seiar0_100[, id,])
matplot(seiar0_100[, t_id,],seiar0_100[, id,], 
        xlab = "Günler", 
        ylab = "Vaka sayısı",
        type = "l", lty = 1, col="grey",
        xlim=c(0,30))
lines(seiar0_100[, 1,1],mean,col="purple")




# Kümülatif NoV ölümleri
id<- which(idx=="cum_deaths_all")
mean <- rowMeans(seiar0_100[, id,])
matplot(seiar0_100[, t_id,],seiar0_100[, id,], 
        xlab = "Günler", 
        ylab = "kümülatif ölüm sayısı",
        type = "l", lty = 1, col="grey",
        xlim=c(0,365))
lines(seiar0_100[, 1,1],mean,col="purple")

# Yaşa göre yığılmış ölümler 

time <- (seiar0_100[, t_id,1])
t<-rep(time,4)
age <-  c(rep(c("0_4") , length(time)), 
          rep(c("5_14") , length(time)),
          rep(c("15_64") , length(time)), 
          rep(c("65+") , length(time)))
deaths <- c( rowMeans(seiar0_100[,which(idx=="cum_deaths[1]") ,]),
             rowMeans(seiar0_100[,which(idx=="cum_deaths[2]") ,]),
             rowMeans(seiar0_100[,which(idx=="cum_deaths[3]") ,]),
             rowMeans(seiar0_100[,which(idx=="cum_deaths[4]") ,]))
df <- data.frame(t,age,deaths)

# Yaşa göre zaman içinde yığılmış kümülatif ölümler
ggplot(df, aes(fill=age, y=deaths, x=t)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("günler")+ ylab("Kümülatif Ölümler")


# Yaşa göre zaman içinde göreceli yığılmış kümülatif ölümler
ggplot(df, aes(fill=age, y=deaths, x=t)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("günler")+ ylab("Kümülatif Ölümlerin Oranı")


#Yaş grubuna göre simüle edilmiş NoV vakaları

# yaşa göre vakalar 
time <- (seiar0_100[, t_id,1])
t<-rep(time,4)
age <-  c(rep(c("0_4") , length(time)), 
          rep(c("5_14") , length(time)),
          rep(c("15_64") , length(time)), 
          rep(c("65+") , length(time)))
cases <- c( 1000*(rowMeans(seiar0_100[,which(idx=="cum_cases[1]") ,])/pop[1]),
            1000*( rowMeans(seiar0_100[,which(idx=="cum_cases[2]") ,])/pop[2]),
            1000*(rowMeans(seiar0_100[,which(idx=="cum_cases[3]") ,])/pop[3]),
            1000*(rowMeans(seiar0_100[,which(idx=="cum_cases[4]") ,]))/pop[4])
df <- data.frame(t,age,cases)

ggplot(df,aes(color=age, y=cases, x=t)) + 
  geom_line() +
  xlab("günler")+ ylab("1000 nüfus başına kümülatif insidans oranı")


# Analiz çıktısı

id<- which(idx=="cum_deaths_all")
base_deaths <- mean(seiar0_100[365 ,id, ]) # bir yıl sonra NoV'den kümülatif ölümler 
print(base_deaths)

id<- which(idx=="cum_cases_all")
base_cases <- mean(seiar0_100[365 ,id, ]) # bir yıl sonra NoV'den kümülatif vakalar 
print(base_cases)

id<- which(idx=="cum_vaccines_all")
base_doses <- mean(seiar0_100[365 ,id, ]) # bir yıl sonra kümülatif aşı dozu (0) 
print(base_doses)


# 3. Yeni bir aşı mevcut hale geliyor -------------------------------------




## Senaryo 1

seiar1 <- seiar_generator$new(
  V_ini = c(1:n_age)*0,
  S_ini = as.numeric(round(N*pop_distr - c(1,0,0,0))),
  E_ini = c(1:n_age)*0,
  I_ini = c(1,0,0,0),
  A_ini = c(1:n_age)*0,
  VA_ini = c(1:n_age)*0,
  R_ini = c(1:n_age)*0,
  N_age = n_age,
  cfr= mu,
  m=transmission,
  rep_ratio=rep_ratio,
  vac_cov= c(0,0,0,0),## <------- Farklı yaş gruplarında aşılama kapsamını seçin (0-1 arası)
  vac_eff = 0,       ## <------- Aşı etkinliğini doldurun (NoVax veya Vomax)
  vac_imm = 0,       ## <------- AB korunmasına geçiş oranını doldurun (NoVax veya Vomax) (1/gün)
  beta = 0           ## <--------2. bölümden en iyi uyumlu betanızla doldurun
)

# Modeli çalıştır
seiar1_100 <- seiar1$run(0:t_end, replicate = 100)

# Analiz çıktısı
id<- which(idx=="cum_deaths_all")
sc1_deaths <- mean(seiar1_100[365 ,id, ]) ## 1 yıl sonra kümülatif ölümler

id<- which(idx=="cum_cases_all")
sc1_cases <- mean(seiar1_100[365 ,id, ]) ## 1 yıl sonra kümülatif vakalar

id<- which(idx=="cum_vaccines_all")
sc1_doses <- mean(seiar1_100[365 ,id, ]) ## 1 yıl sonra kümülatif doz

sc1_avdeaths<-base_deaths-sc1_deaths ## 1 yıl sonra önlenen ölümler

sc1_avcases<-base_cases-sc1_cases    ## 1 yıl sonra önlenen vakalar

dose_per_avdeath1<- sc1_doses/sc1_avdeaths ## 1 ölümü önlemek için gereken doz sayısı 

dose_per_avcase1<- sc1_doses/sc1_avcases ## 1 vakayı önlemek için gereken doz sayısı 


print(paste("1 yıl sonra kümülatif ölümler",sc1_deaths))
print(paste("1 yıl sonra kümülatif vakalar",sc1_cases))
print(paste("1 yıl sonra kümülatif doz",sc1_doses))

print(paste("1 yıl sonra önlenen ölümler",sc1_avdeaths))
print(paste("1 yıl sonra önlenen vakalar",sc1_avcases))
print(paste("1 ölümü önlemek için gereken doz sayısı",dose_per_avdeath1))
print(paste("1 vakayı önlemek için gereken doz sayısı",dose_per_avcase1))



## Senaryo 2

seiar2 <- seiar_generator$new(
  V_ini = c(1:n_age)*0,
  S_ini = as.numeric(round(N*pop_distr - c(1,0,0,0))),
  E_ini = c(1:n_age)*0,
  I_ini = c(1,0,0,0),
  A_ini = c(1:n_age)*0,
  VA_ini = c(1:n_age)*0,
  R_ini = c(1:n_age)*0,
  N_age = n_age,
  cfr= mu,
  m=transmission,
  rep_ratio=rep_ratio,
  vac_cov= c(0,0,0,0),## <------- Farklı yaş gruplarında aşılama kapsamını seçin (0-1 arası)
  vac_eff = 0,       ## <------- Aşı etkinliğini doldurun (NoVax veya Vomax)
  vac_imm = 0,       ## <------- AB korunmasına geçiş oranını doldurun (NoVax veya Vomax) (1/gün)
  beta = 0           ## <--------2. bölümden en iyi uyumlu betanızla doldurun
)

# Modeli çalıştır
seiar2_100 <- seiar2$run(0:t_end, replicate = 100)

# Analiz çıktısı
id<- which(idx=="cum_deaths_all")
sc2_deaths <- mean(seiar2_100[365 ,id, ]) ## 1 yıl sonra kümülatif ölümler

id<- which(idx=="cum_cases_all")
sc2_cases <- mean(seiar2_100[365 ,id, ]) ## 1 yıl sonra kümülatif vakalar

id<- which(idx=="cum_vaccines_all")
sc2_doses <- mean(seiar2_100[365 ,id, ]) ## 1 yıl sonra kümülatif doz

sc2_avdeaths<-base_deaths-sc2_deaths ## 1 yıl sonra önlenen ölümler

sc2_avcases<-base_cases-sc2_cases    ## 1 yıl sonra önlenen vakalar

dose_per_avdeath2<- sc2_doses/sc2_avdeaths ## 1 ölümü önlemek için gereken doz sayısı 

dose_per_avcase2<- sc2_doses/sc2_avcases ## 1 vakayı önlemek için gereken doz sayısı 

print(paste("1 yıl sonra kümülatif ölümler",sc2_deaths))
print(paste("1 yıl sonra kümülatif vakalar",sc2_cases))
print(paste("1 yıl sonra kümülatif doz",sc2_doses))

print(paste("1 yıl sonra önlenen ölümler",sc2_avdeaths))
print(paste("1 yıl sonra önlenen vakalar",sc2_avcases))
print(paste("1 ölümü önlemek için gereken doz sayısı",dose_per_avdeath2))
print(paste("1 vakayı önlemek için gereken doz sayısı",dose_per_avcase2))



## Senaryo 3

seiar3 <- seiar_generator$new(
  V_ini = c(1:n_age)*0,
  S_ini = as.numeric(round(N*pop_distr - c(1,0,0,0))),
  E_ini = c(1:n_age)*0,
  I_ini = c(1,0,0,0),
  A_ini = c(1:n_age)*0,
  VA_ini = c(1:n_age)*0,
  R_ini = c(1:n_age)*0,
  N_age = n_age,
  cfr= mu,
  m=transmission,
  rep_ratio=rep_ratio,
  vac_cov= c(0,0,0,0),## <------- Farklı yaş gruplarında aşılama kapsamını seçin (0-1 arası)
  vac_eff = 0,       ## <------- Aşı etkinliğini doldurun (NoVax veya Vomax)
  vac_imm = 0,       ## <------- AB korunmasına geçiş oranını doldurun (NoVax veya Vomax) (1/gün)
  beta = 0           ## <--------2. bölümden en iyi uyumlu betanızla doldurun
)

# Modeli çalıştır
seiar3_100 <- seiar3$run(0:t_end, replicate = 100)

# Analiz çıktısı
id<- which(idx=="cum_deaths_all")
sc3_deaths <- mean(seiar3_100[365 ,id, ]) ## 1 yıl sonra kümülatif ölümler

id<- which(idx=="cum_cases_all")
sc3_cases <- mean(seiar3_100[365 ,id, ]) ## 1 yıl sonra kümülatif vakalar

id<- which(idx=="cum_vaccines_all")
sc3_doses <- mean(seiar3_100[365 ,id, ]) ## 1 yıl sonra kümülatif doz

sc3_avdeaths<-base_deaths-sc3_deaths ## 1 yıl sonra önlenen ölümler

sc3_avcases<-base_cases-sc3_cases    ## 1 yıl sonra önlenen vakalar

dose_per_avdeath3<- sc3_doses/sc3_avdeaths ## 1 ölümü önlemek için gereken doz sayısı 

dose_per_avcase3<- sc3_doses/sc3_avcases ## 1 vakayı önlemek için gereken doz sayısı 

print(paste("1 yıl sonra kümülatif ölümler",sc3_deaths))
print(paste("1 yıl sonra kümülatif vakalar",sc3_cases))
print(paste("1 yıl sonra kümülatif doz",sc3_doses))

print(paste("1 yıl sonra önlenen ölümler",sc3_avdeaths))
print(paste("1 yıl sonra önlenen vakalar",sc3_avcases))
print(paste("1 ölümü önlemek için gereken doz sayısı",dose_per_avdeath3))
print(paste("1 vakayı önlemek için gereken doz sayısı",dose_per_avcase3))

