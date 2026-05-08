library(multcompView)
df<-read.csv('data_meanGrowthDifference.csv')
df$logVals<-sign(df$Growth.difference) * log(1 + abs(df$Growth.difference))
df$HostGroup<-as.factor(df$HostGroup)
df$PhyloGroup<-as.factor(df$PhyloGroup)

L<-lm('logVals~HostGroup',data=df)
print(shapiro.test(residuals(L)))
A<-aov(L)
print(summary(A))
print('#######################')


L<-lm('logVals~PhyloGroup',data=df)
print(shapiro.test(residuals(L)))
A<-aov(L)
print(summary(A))
T<-TukeyHSD(A)
print(T)
print(multcompLetters4(L,T))

L<-lm('logVals~Continent.of.isolation',data=df)
print(shapiro.test(residuals(L)))
A<-aov(L)
print(summary(A))
T<-TukeyHSD(A)
print(T)
print(multcompLetters4(L,T))

