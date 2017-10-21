library(R.matlab)
library(dplyr)
library(ggplot2)

data <- readMat("olympics.mat")
str(data)
(male100 <- data$male100[ ,1:2])
(female100 <- data$female100[ ,1:2])
(male200 <- data$male200[ ,1:2])
(female200 <- data$female200[ ,1:2])
(male400 <- data$male400[ ,1:2])
(female400 <- data$female400[ ,1:2])

# 2012 men's 400    43.94
# 2016 men's 400    43.03
# 2012 women's 400  49.55
# 2016 women's 400  49.44

# 2012 men's 200    19.32
# 2016 men's 200    19.78
# 2012 women's 200  21.88
# 2016 women's 200  21.78

# 2012 men's 100     9.63
# 2016 men's 100     9.81
# 2012 women's 100  10.75
# 2016 women's 100  10.71


### men's 400
male400 <- as.data.frame(male400) %>% rename(year = V1, seconds = V2)
male400 <- male400 %>% bind_rows(
  c(year = 2012, seconds = 43.94), c(year = 2016, seconds = 43.03))
ggplot(male400, aes(x = year, y = seconds)) + geom_line() + ggtitle("Men's 400m Olympic winning times 1896-2016")

male400_1996 <- male400 %>% filter(year < 1997)
ggplot(male400_1996, aes(x = year, y = seconds)) + geom_line() + ggtitle("Men's 400m Olympic winning times 1896-1996")

### women's 400
female400 <- as.data.frame(female400) %>% rename(year = V1, seconds = V2)
female400 <- female400 %>% bind_rows(
  c(year = 2012, seconds = 49.55), c(year = 2016, seconds = 49.44))
ggplot(female400, aes(x = year, y = seconds)) + geom_line() + ggtitle("Women's 400m Olympic winning times 1964-2016")

female400_1996 <- female400 %>% filter(year < 1997)
ggplot(female400_1996, aes(x = year, y = seconds)) + geom_line() + ggtitle("Women's 400m Olympic winning times 1964-1996")


### men's 200
male200 <- as.data.frame(male200) %>% rename(year = V1, seconds = V2)
male200 <- male200 %>% bind_rows(
  c(year = 2012, seconds = 19.32), c(year = 2016, seconds = 19.78))
ggplot(male200, aes(x = year, y = seconds)) + geom_line() + ggtitle("Men's 200m Olympic winning times 1896-2016")

male200_1996 <- male200 %>% filter(year < 1997)
ggplot(male200_1996, aes(x = year, y = seconds)) + geom_line() + ggtitle("Men's 200m Olympic winning times 1896-1996")

### women's 200
female200 <- as.data.frame(female200) %>% rename(year = V1, seconds = V2)
female200 <- female200 %>% bind_rows(
  c(year = 2012, seconds = 21.88), c(year = 2016, seconds = 21.78))
ggplot(female200, aes(x = year, y = seconds)) + geom_line() + ggtitle("Women's 200m Olympic winning times 1948-2016")

female200_1996 <- female200 %>% filter(year < 1997)
ggplot(female200_1996, aes(x = year, y = seconds)) + geom_line() + ggtitle("Women's 200m Olympic winning times 1948-1996")


### men's 100
male100 <- as.data.frame(male100) %>% rename(year = V1, seconds = V2)
male100 <- male100 %>% bind_rows(
  c(year = 2012, seconds = 9.63), c(year = 2016, seconds = 9.81))
ggplot(male100, aes(x = year, y = seconds)) + geom_line() + ggtitle("Men's 100m Olympic winning times 1896-2016")

male100_1996 <- male100 %>% filter(year < 1997)
ggplot(male100_1996, aes(x = year, y = seconds)) + geom_line() + ggtitle("Men's 100m Olympic winning times 1896-1996")

### women's 100
female100 <- as.data.frame(female100) %>% rename(year = V1, seconds = V2)
female100 <- female100 %>% bind_rows(
  c(year = 2012, seconds = 10.75), c(year = 2016, seconds = 10.71))
ggplot(female100, aes(x = year, y = seconds)) + geom_line() + ggtitle("Women's 100m Olympic winning times 1928-2016")

female100_1996 <- female100 %>% filter(year < 1997)
ggplot(female100_1996, aes(x = year, y = seconds)) + geom_line() + ggtitle("Women's 100m Olympic winning times 1928-1996")


