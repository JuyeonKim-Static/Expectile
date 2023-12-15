Spatial Quantile Correlation
================
JuyeonKim
2023-09-19

``` r
expec_corr_seoul <- read_csv("expec_corr_seoul.csv")
```

    ## New names:
    ## Rows: 15 Columns: 10
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (9): 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
expec_corr_seoul_gangwon <- read_csv("expec_corr_Gangwon.csv")
```

    ## New names:
    ## Rows: 5 Columns: 10
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (9): 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
expec_corr_seoul <- rbind(expec_corr_seoul, expec_corr_seoul_gangwon[c(1, 4, 5),])
region_coord <- read_excel("region_coord.xlsx")
region_coord_gangwon <- read_excel("gangwonregion_coord.xlsx")
expec_corr_seoul[c('lat', 'long')] <- NA
expec_corr_seoul[c('lat', 'long')][1:15,] <- region_coord[c('Latitude', 'Longitude')][1:15,]
expec_corr_seoul[c('lat', 'long')][16:18,] <- region_coord_gangwon[c('Latitude', 'Longitude')][c(4, 5, 3),]
names(expec_corr_seoul)[1] <- 'Loc'
expec_corr_seoul
```

    ## # A tibble: 18 × 12
    ##    Loc         `0.1` `0.2` `0.3` `0.4` `0.5` `0.6` `0.7` `0.8` `0.9`   lat  long
    ##    <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1 Gangwon     0.548 0.611 0.659 0.701 0.741 0.781 0.826 0.884 0.978  37.7  129.
    ##  2 Seoul       1     1     1     1     1     1     1     1     1      37.6  127.
    ##  3 Gyeonggi    0.811 0.852 0.883 0.909 0.933 0.957 0.985 1.02  1.08   37.3  127.
    ##  4 Chungcheon… 0.500 0.571 0.622 0.666 0.709 0.754 0.804 0.869 0.971  36.2  128.
    ##  5 Gyeongsang… 0.533 0.604 0.655 0.697 0.739 0.784 0.836 0.900 1.00   36.6  129.
    ##  6 Daegu       0.464 0.534 0.584 0.628 0.670 0.713 0.762 0.827 0.934  35.9  129.
    ##  7 Jeollabuk   0.637 0.691 0.732 0.768 0.802 0.837 0.877 0.927 1.00   35.8  127.
    ##  8 Ulsan       0.404 0.468 0.515 0.557 0.600 0.645 0.697 0.765 0.880  35.6  129.
    ##  9 Gwangju     0.513 0.582 0.633 0.677 0.719 0.763 0.811 0.871 0.969  35.2  127.
    ## 10 Busan       0.375 0.439 0.489 0.535 0.579 0.626 0.683 0.759 0.890  35.1  129.
    ## 11 Jeollanam   0.405 0.476 0.527 0.572 0.615 0.659 0.706 0.769 0.872  34.5  126.
    ## 12 Jeju        0.325 0.385 0.433 0.477 0.523 0.576 0.642 0.734 0.894  33.3  126.
    ## 13 Gyeongsang… 0.406 0.481 0.535 0.583 0.629 0.679 0.739 0.820 0.958  35.2  128.
    ## 14 Incheon     0.742 0.795 0.829 0.859 0.886 0.913 0.944 0.982 1.04   37.7  126.
    ## 15 Chungcheon… 0.699 0.751 0.789 0.820 0.849 0.878 0.910 0.949 1.01   36.8  127.
    ## 16 Sockcho     0.475 0.543 0.593 0.638 0.680 0.724 0.773 0.835 0.927  38.2  129.
    ## 17 Chuncheon   0.726 0.773 0.805 0.833 0.858 0.882 0.907 0.936 0.974  37.9  128.
    ## 18 Yeongwol    0.560 0.628 0.678 0.720 0.761 0.803 0.850 0.911 1.01   37.2  128.

``` r
quant_corr_seoul <- read_csv("quant_corr_seoul.csv")
```

    ## New names:
    ## Rows: 15 Columns: 10
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (9): 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
quant_corr_seoul_gangwon <- read_csv("quant_corr_Gangwon.csv")
```

    ## New names:
    ## Rows: 5 Columns: 10
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (9): 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
quant_corr_seoul <- rbind(quant_corr_seoul, quant_corr_seoul_gangwon[c(1, 4, 5),])
region_coord <- read_excel("region_coord.xlsx")
region_coord_gangwon <- read_excel("gangwonregion_coord.xlsx")
quant_corr_seoul[c('lat', 'long')] <- NA
quant_corr_seoul[c('lat', 'long')][1:15,] <- region_coord[c('Latitude', 'Longitude')][1:15,]
quant_corr_seoul[c('lat', 'long')][16:18,] <- region_coord_gangwon[c('Latitude', 'Longitude')][c(4, 5, 3),]
names(quant_corr_seoul)[1] <- 'Loc'
quant_corr_seoul
```

    ## # A tibble: 18 × 12
    ##    Loc         `0.1` `0.2` `0.3` `0.4` `0.5` `0.6` `0.7` `0.8` `0.9`   lat  long
    ##    <chr>       <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>
    ##  1 Gangwon     0.496 0.551 0.602 0.643 0.691 0.737 0.804 0.845 0.931  37.7  129.
    ##  2 Seoul       1     1     1     1     1     1     1     1     1      37.6  127.
    ##  3 Gyeonggi    0.794 0.842 0.877 0.903 0.922 0.948 0.983 1.01  1.05   37.3  127.
    ##  4 Chungcheon… 0.450 0.519 0.585 0.623 0.667 0.720 0.774 0.849 0.974  36.2  128.
    ##  5 Gyeongsang… 0.532 0.582 0.628 0.671 0.712 0.762 0.823 0.904 0.992  36.6  129.
    ##  6 Daegu       0.463 0.522 0.581 0.633 0.674 0.714 0.754 0.814 0.944  35.9  129.
    ##  7 Jeollabuk   0.579 0.638 0.691 0.730 0.772 0.815 0.858 0.904 1.00   35.8  127.
    ##  8 Ulsan       0.398 0.451 0.499 0.537 0.588 0.644 0.691 0.775 0.860  35.6  129.
    ##  9 Gwangju     0.472 0.544 0.601 0.657 0.695 0.737 0.798 0.868 0.959  35.2  127.
    ## 10 Busan       0.350 0.385 0.431 0.487 0.528 0.586 0.633 0.722 0.822  35.1  129.
    ## 11 Jeollanam   0.378 0.417 0.462 0.505 0.558 0.611 0.654 0.725 0.824  34.5  126.
    ## 12 Jeju        0.265 0.326 0.356 0.388 0.416 0.477 0.530 0.584 0.717  33.3  126.
    ## 13 Gyeongsang… 0.433 0.482 0.524 0.567 0.613 0.658 0.718 0.801 0.900  35.2  128.
    ## 14 Incheon     0.693 0.753 0.798 0.834 0.867 0.892 0.930 0.982 1.05   37.7  126.
    ## 15 Chungcheon… 0.663 0.719 0.770 0.800 0.830 0.865 0.904 0.948 1.02   36.8  127.
    ## 16 Sockcho     0.438 0.506 0.535 0.583 0.636 0.690 0.739 0.813 0.938  38.2  129.
    ## 17 Chuncheon   0.671 0.724 0.746 0.779 0.820 0.850 0.895 0.933 0.971  37.9  128.
    ## 18 Yeongwol    0.525 0.592 0.640 0.687 0.729 0.777 0.828 0.903 1.01   37.2  128.

``` r
x <- cbind(quant_corr_seoul$long, 
           quant_corr_seoul$lat)
y <- quant_corr_seoul$`0.1`

mKrigObject<- mKrig( x,y,lambda=.00,cov.args= list(     aRange= 5.87,
                                    Covariance="Matern",
                                    smoothness=1.0),
                                        sigma2=.157
                                 )
gridList<- list( x = seq(124, 130, length.out = 100),
                 y = seq( 33,   39, length.out = 100)) 

gHat<- predictSurface( mKrigObject, gridList)
 
 image.plot(gHat$x, gHat$y,  gHat$z  ,zlim =c(0.4,1.1), xlab = "Longitude", ylab = "Latitude" )
 world( add=TRUE )
 title("0.1 Quantile Correlation Coefficients")
```

![patialQCorr_files/figure-gfm/unnamed-chunk-3-1.png](https://github.com/JuyeonKim-Static/Expectile/blob/d190562688b62e429a8281d8a2dd97b456b963fe/PM10_data/image/unnamed-chunk-3-1.png))<!-- -->
``` r
x <- cbind(quant_corr_seoul$long, 
           quant_corr_seoul$lat)
y <- quant_corr_seoul$`0.9`

mKrigObject<- mKrig( x,y,
                    lambda=.00,
                    cov.args= list(     aRange= 5.87,
                                    Covariance="Matern",
                                    smoothness=1.0),
                                        sigma2=.157
                                 )
gridList<- list( x = seq(124, 130, length.out = 100),
                 y = seq( 33,   39, length.out = 100)) 

gHat<- predictSurface( mKrigObject, gridList)
 
 image.plot(gHat$x, gHat$y,  gHat$z, zlim =c(0.4,1.1), xlab = "Longitude", ylab = "Latitude" )
 world( add=TRUE )
 title("0.9 Quantile Correlation Coefficients")
```
![patialQCorr_files/figure-gfm/unnamed-chunk-3-2.png](https://github.com/JuyeonKim-Static/Expectile/blob/d190562688b62e429a8281d8a2dd97b456b963fe/PM10_data/image/unnamed-chunk-3-2.png))<!-- -->

``` r
x <- cbind(expec_corr_seoul$long, 
           expec_corr_seoul$lat)
y <- expec_corr_seoul$`0.1`

mKrigObject<- mKrig( x,y,
                    lambda=.00,
                    cov.args= list(     aRange= 5.87,
                                    Covariance="Matern",
                                    smoothness=1.0),
                                        sigma2=.157
                                 )
gridList<- list( x = seq(124, 130, length.out = 100),
                 y = seq( 33,   39, length.out = 100)) 

gHat<- predictSurface( mKrigObject, gridList)
 
 image.plot(gHat$x, gHat$y,  gHat$z  ,zlim =c(0.4,1.1), xlab = "Longitude", ylab = "Latitude" )
 world( add=TRUE )
 title("0.1 Expectile Correlation Coefficients")
```
![patialQCorr_files/figure-gfm/unnamed-chunk-4-1.png](https://github.com/JuyeonKim-Static/Expectile/blob/d190562688b62e429a8281d8a2dd97b456b963fe/PM10_data/image/unnamed-chunk-4-1.png))<!-- -->

``` r
x <- cbind(expec_corr_seoul$long, 
           expec_corr_seoul$lat)
y <- expec_corr_seoul$`0.9`

mKrigObject<- mKrig( x,y,
                    lambda=.00,
                    cov.args= list(     aRange= 5.87,
                                    Covariance="Matern",
                                    smoothness=1.0),
                                        sigma2=.157
                                 )
gridList<- list( x = seq(124, 130, length.out = 100),
                 y = seq( 33,   39, length.out = 100)) 

gHat<- predictSurface( mKrigObject, gridList)
 
 image.plot(gHat$x, gHat$y,  gHat$z, zlim =c(0.4,1.1), xlab = "Longitude", ylab = "Latitude" )
 world( add=TRUE )
 title("0.9 Expectile Correlation Coefficients")
```

![patialQCorr_files/figure-gfm/unnamed-chunk-4-2.png](https://github.com/JuyeonKim-Static/Expectile/blob/d190562688b62e429a8281d8a2dd97b456b963fe/PM10_data/image/unnamed-chunk-4-2.png))<!-- -->

``` r
library('latex2exp')
library(geosphere)
dist_seoul <- expec_corr_seoul %>% select(c(long, lat)) %>% 
  distm(expec_corr_seoul %>% filter(Loc == 'Seoul') %>% select(c(long, lat)), fun = distHaversine)
dist_seoul <- dist_seoul / 1000 # make distance as km 
plot(dist_seoul, expec_corr_seoul$`0.1`, pch = 16, col = "blue", 
     xlab = "Distance (km)", ylab = TeX("$\\rho_{\\tau}^{(e)}$") ,ylim=c(0.3,1.1),cex.lab= 0.7)
points(dist_seoul, expec_corr_seoul$`0.9`, pch = 16, col = "red")
legend("bottomleft", legend = c(TeX("$\\tau = 0.1$"), TeX("$\\tau = 0.9$")), pch = 16, col = c("blue", "red"))
title(expression(paste(rho[tau]^{(e)}, " with Distance from Seoul")))
```

![patialQCorr_files/figure-gfm/unnamed-chunk-5-1.png](https://github.com/JuyeonKim-Static/Expectile/blob/d190562688b62e429a8281d8a2dd97b456b963fe/PM10_data/image/unnamed-chunk-5-1.png))<!-- -->

``` r
library(geosphere)
dist_seoul <- quant_corr_seoul %>% select(c(long, lat)) %>% 
  distm(quant_corr_seoul %>% filter(Loc == 'Seoul') %>% select(c(long, lat)), fun = distHaversine)
dist_seoul <- dist_seoul / 1000 # make distance as km 
plot(dist_seoul, quant_corr_seoul$`0.1`, pch = 16, col = "blue", 
     xlab = "Distance (km)", ylab = TeX("$\\rho_{\\tau}^{(q)}$") ,ylim=c(0.3,1.1), cex.lab= 0.7)
points(dist_seoul, quant_corr_seoul$`0.9`, pch = 16, col = "red")
legend("bottomleft", legend = c(TeX("$\\tau = 0.1$"), TeX("$\\tau = 0.9$")) , pch = 16, col = c("blue", "red"))
title(expression(paste(rho[tau]^{(q)}, " with Distance from Seoul")))
```

![patialQCorr_files/figure-gfm/unnamed-chunk-6-1.png](https://github.com/JuyeonKim-Static/Expectile/blob/d190562688b62e429a8281d8a2dd97b456b963fe/PM10_data/image/unnamed-chunk-6-1.png))<!-- -->
