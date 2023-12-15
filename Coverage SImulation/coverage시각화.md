coverage_시각화
================
Juyeon Kim
2023-12-15

``` r
library(readxl)
library(ggplot2)
#n=500_rho=0.3
df<- read_excel("신뢰구간길이_new.xlsx")
df2 <- read_excel("coverage_시각화.xlsx")
a<- rep(1:24)
#plot
p <- ggplot(df, aes(x = 1, y = quant_ci1)) +
  geom_violin(fill = "red") +
  labs(title = "Violinplot of quant_ci and expec_ci") +
  ylim(0, 1) +
  geom_violin(aes(x = 2, y = expec_ci1), fill = "blue") +
  geom_violin(aes(x = 3, y = quant_ci2), fill = "red") +
  geom_violin(aes(x = 4, y = expec_ci2), fill = "blue") +
  geom_violin(aes(x = 5, y = quant_ci3), fill = "red") +
  geom_violin(aes(x = 6, y = expec_ci3), fill = "blue") +
  geom_violin(aes(x = 7, y = quant_ci4), fill = "red") +
  geom_violin(aes(x = 8, y = expec_ci4), fill = "blue") +
  geom_violin(aes(x = 9, y = quant_ci5), fill = "red") +
  geom_violin(aes(x = 10, y = expec_ci5), fill = "blue") +
  geom_violin(aes(x = 11, y = quant_ci6), fill = "red") +
  geom_violin(aes(x = 12, y = expec_ci6), fill = "blue") +
  geom_violin(aes(x = 13, y = quant_ci7), fill = "red") +
  geom_violin(aes(x = 14, y = expec_ci7), fill = "blue") +
  geom_violin(aes(x = 15, y = quant_ci8), fill = "red") +
  geom_violin(aes(x = 16, y = expec_ci8), fill = "blue") +
  geom_violin(aes(x = 17, y = quant_ci9), fill = "red") +
  geom_violin(aes(x = 18, y = expec_ci9), fill = "blue") +
  geom_violin(aes(x = 19, y = quant_ci10), fill = "red") +
  geom_violin(aes(x = 20, y = expec_ci10), fill = "blue") +
  geom_violin(aes(x = 21, y = quant_ci11), fill = "red") +
  geom_violin(aes(x = 22, y = expec_ci11), fill = "blue") +
  geom_violin(aes(x = 23, y = quant_ci12), fill = "red") +
  geom_violin(aes(x = 24, y = expec_ci12), fill = "blue") +
  # geom_violin(aes(x = 25, y = quant_ci13), fill = "red") +
  # geom_violin(aes(x = 26, y = expec_ci13), fill = "blue") +
  # geom_violin(aes(x = 27, y = quant_ci14), fill = "red") +
  # geom_violin(aes(x = 28, y = expec_ci14), fill = "blue") +
  # geom_violin(aes(x = 29, y = quant_ci15), fill = "red") +
  # geom_violin(aes(x = 30, y = expec_ci15), fill = "blue") +
  geom_vline(xintercept=c(6.5, 12.5, 18.5), linetype = "dashed", col = "gray30")+
  scale_x_continuous(breaks = c(3.5, 9.5, 15.5, 21.5), 
                     labels = c("df = 3", "df = 5", "df = 7", "MVN"))
p <- p + geom_point(data = df2[1:24,], aes(x = a[1:24], y = coverage_mean[1:24]), color = "green3", size = 2 ,shape=3)
# Adjust the y-axis labels and limits for the scatter plot
p <- p + scale_y_continuous(sec.axis = sec_axis(~., name = "coverage_mean", breaks = c(0.85, 0.899, 0.856))) +ylim(0,1)
```

    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.
    ## Scale for y is already present.
    ## Adding another scale for y, which will replace the existing scale.

``` r
p <- p + xlab("degree of Freedom") + ylab("Length of CI") +theme_bw()

p
```

    ## Warning: Removed 3 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 20 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 4 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 5 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 25 rows containing non-finite values (`stat_ydensity()`).

    ## Warning: Removed 1 rows containing non-finite values (`stat_ydensity()`).
    ## Removed 1 rows containing non-finite values (`stat_ydensity()`).

![](https://github.com/JuyeonKim-Static/Expectile/blob/b39d4a079ecf22d1356153c0915a0680e9d19c12/Coverage%20SImulation/violinplot.png)<!-- -->
