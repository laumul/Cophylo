main <- function(){
names <- read.table("names.txt")

for (i in 1:length(names$V1)){
  
fl <- paste0("brn/", names[i,])
data <- read.table(fl)
ln <- length(data$V1) + 1

a <- mean(data[["V3"]])
b <- mean(data[["V4"]])
c <- mean(data[["V5"]])
d <- mean(data[["V6"]])
e <- mean(data[["V7"]])

data[ln,]$V1 <- "S_terminal_branches"
data[ln,]$V2 <- data[2,]$V2
data[ln,]$V3 <- a
data[ln,]$V4 <- b
data[ln,]$V5 <- c
data[ln,]$V6 <- d
data[ln,]$V7 <- e
nm <- paste0("aver/", names[i,])
write.table(data, nm )

}
}
main ()