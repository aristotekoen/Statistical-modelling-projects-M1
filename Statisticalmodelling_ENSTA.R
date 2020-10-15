

#Test de Shapiro
shapiro.test(data$Debit)
?Cauchy
#Q1: Calcul de P(X>1):

1-pcauchy(1,0,1)
#0.25 La loi de Cauchy est plus étalée que la loi normale, donc plus de poids sur evenements rares

1- pnorm(1,0,1)
#0.1586553

#plot de la fonction de répartition d'une Cauchy(0,1) vs loi Normale centrée réduite
curve(pcauchy(x,0,1),from = c(-3,3))
curve(pnorm(x,0,1),add = TRUE)
?curve
?cumsum

sample = rcauchy(n = 500)
res = cumsum(sample)/seq(1,500) 
#chaque composante du vecteur est la moyenne empirique sur les observations précédentes.
#On observe un comportement aléatoire et n'observe à priori aucune convergence de la moyenne empirique vers une certaine valeur. 
#En effet la loi des grands nombres n'est pas applicable sachant que son esperance n'est pas définie. 
par(mfrow=c(1,1))
plot(seq(1,500),res,type = 'l')

#Partie 2:
data = read.csv('abermule.csv')

#On initialise le germe à 1
set.seed(1)

#On génère un echantillon de loi de Cauchy de paramètre 0,1 de taille 3
ech1=rcauchy(3,0,1)

#On crée une fonction qui calcule la vraisemblance d'une Cauchy 0,1
log_a= function(x,echant)
{
  n=length(echant)
  y=-n*log(pi) - sum(log(1+(echant-x)^2))
  return(y)
}

int1=seq(-5,5,0.01)
int2=seq(-500,500,1)
#plot de la fonction de vraisemblance sur une fenêtre -5,5. Il y a un max local au voisinage de -4
par(mfrow=c(1,1))

plot(int1,lapply(int1,log_a,echant=ech1),type='l',col='blue',xlab = 'a',ylab = 'log_l(a)',main = 'Log Vraisemblance sur [-5,5]')
legend(x=-5,y=-8,cex=0.8,legend=c('log vraisemblance'),lty=1,col='blue')

plot(int2,lapply(int2,log_a,echant=ech1),type='l',col='blue',xlab='a',ylab='log_l(a)',main='Log Vraisemblance sur [-500,500]')
legend(x=-500,y=-8,cex=0.8,legend=c('log vraisemblance'),lty=1,col='blue')

abline(v=1)
#Il suffit de choisir un point d'initialisation positif et grand afin de s'assurer de ne pas atteindre le max local 

#Question 2:

#nlm fait appel à un algorithme d'optimisation de type newton.
IQR(data$Debit)
B=50
ech50=replicate(B,rcauchy(100,0,1),simplify = TRUE)


#on crée la fonction calculant la fonction de vraisemblance en fonction de a et b
log_a_b= function(x,echant){
   a=x[1]; b=x[2]
    n=length(echant)
    res=n*log(pi*b) + sum(log(1+((echant-a)/b)^2))
    
    attr(res,"gradient")=c(-2*sum((echant-a)/((b)^2 + (echant-a)^2)),(n/b)-2*sum((echant-a)^2/((b^3)+((echant-a)^2)*b)))
    
     res}

#création d'une matrice avec les valeurs de la hessienne en le maximum de vraisemblance
matrice=matrix(ncol=2,nrow = 50)
for(i in seq(1,50)){
  c=nlm(log_a_b,c(50,50),echant=ech50[,i])["estimate"]
  matrice[i,1]=c[[1]][1]
  matrice[i,2]=c[[1]][2]
  
}

#plot des histogrammes des EMV
par(mfrow=c(1,1))
hist(matrice[,1],breaks=6,proba=TRUE, ylab = 'nb d échantillons',xlab = 'EMV de a',main=paste("Histogramme des EMV de a"))
abline(v=mean(matrice[,1]),col='red')
legend(0.12,3,legend=c(paste('moyenne des EMV de a=',(round(mean(matrice[,1]),2))),paste('variance des EMV de a=',(round(var(matrice[,1]),2)))),lty=1,col=c('red','blue'),cex=0.5)
hist(matrice[,2],proba=TRUE,breaks=7, ylab = 'nb d échantillons',xlab = 'EMV de b',main=paste("Histogramme des EMV de b"))
abline(v=mean(matrice[,2]),col='red')
legend(1.01,3,c(paste('moyenne des EMV de b=',(round(mean(matrice[,2]),2))),paste('variance des EMV de b=',(round(var(matrice[,2]),2)))),lty=1,col=c('red','blue'),cex=0.5)


#plot de l'histogramme du jeu de données
rain= c(nlm(log_a_b,c(50,50),data$Debit,hessian = TRUE)['estimate'][[1]][1],nlm(log_a_b,c(50,50),data$Debit,hessian = TRUE)['estimate'][[1]][2])
hist(data$Debit,breaks=20,proba=TRUE,xlab='Débit', main='Histogramme du débit')
curve(dcauchy(x,rain[1],rain[2]),from = -100, to=1000,add=TRUE,col='red')

#On observe que la loi de Cauchy du paramètre estimé correspond bien àl'histogramme de notre jeu de données. De plus la médiane de l'échantillon est de 191.679 ce qui est très proche de l'estimation de notre paramètre a par la méthode du maximum de vraisemblance qui est égale à 190.1 
median(data$Debit)

#supponson que le modèle est Gaussien, alors l'estimateur du max de vraisemblance est la moyenne empirique et la variance empirique biaisée
taille=length(data$Debit)
m=mean(data$Debit)
v=var(data$Debit)*(taille-1)/taille

#superpositiion de la loi normale deparamètre l'estimateur du max de vraisemblance pour le modèle Gaussien
curve(dnorm(x,m,sqrt(v)),add=TRUE,col='blue')
legend(400,0.012,c('Cauchy EMV','Gaussienne EMV'),cex = 0.6,col=c('red','blue'),lty = 1)
boxplot(data$Debit)


#La densité de la loi normale pour les estimateurs du max de vraisemblance du jeu de données est bien plus aplatie et donne donc un poids plus faible aux valeurs de haut débit
#la loi de Cauchy étant bien plus ressérée elle donne un poids plus important aux valeurs autour de la mediane de l'échantillon.  
legend()
hess=nlm(log_a_b,c(50,50),echant=data$Debit,hessian = TRUE)['hessian']
hessian= matrix(ncol=2,nrow=2)
hessian[1,1]=hess[[1]][1,1]
hessian[1,2]=hess[[1]][1,2]
hessian[2,1]=hess[[1]][2,1]
hessian[2,2]=hess[[1]][2,2]


#On calcule l'inverse de la hessienne pour obtenir l'estimation de notre matrice de variance

solve(hessian)

#Calcul de la p valeur
2-2*(pnorm(sqrt(hessian[2,2])*(rain[2]-1))) #la p valeur est quasi nulle donc on rejette H0 ce qui justifie l'importance de prendre un paramètre d'echelle plus important


curve(dchisq(x,df=1),from=-0, to=10)
curve(dcauchy(x,190,1),from=0,to =600)
cst=sqrt(hessian[2,2])*(rain[2]-1)

proba=1:(length(data$Debit))/(length(data$Debit)+1)
sorted=sort(data$Debit)
plot(qcauchy(proba,rain[1],rain[2]),sorted)
