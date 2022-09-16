dyn.load("itersolve")

n <- 100
a <- 3
A <- matrix(0,n,n)
A[row(A) == (col(A) + 1)] <- - 1
A[col(A) == (row(A) + 1)] <- - 1
diag(A) <- a 

v <- rep(c(1,0), n/2)           
b <- A%*%v
x0 <- rep(0,n)

t <- rep(0,100)
it <- rep(0,100)
error <- rep(0,100)

system.time(a1 <- .C("jacobi", 
		     as.double(A), 
		     as.double(b), 
		     x0 = as.double(x0), 
		     it = as.double(it), 
		     t = as.double(t),
		     error = as.double(error),
		     as.double(v)))

system.time(a2 <- .C("jacobi_parallel", 
		     as.double(A), 
		     as.double(b), 
		     x0 = as.double(x0), 
		     it = as.double(it), 
		     t = as.double(t),
		     error = as.double(error),
		     as.double(v)))


system.time(a3 <- .C("gauss_seidel", 
		     as.double(A), 
		     as.double(b), 
		     x0 = as.double(x0), 
		     it = as.double(it), 
		     t = as.double(t),
		     error = as.double(error),
		     as.double(v)))




pdf(file="plot_iter_2.pdf")
plot(a1$it, a1$error, type = "l", xlab = "Iterations", ylab = "Relative Error",  col = "red")
lines(a2$it, a2$error, type = "l", xlab = "Iterations", ylab = "Relative Error",  col = "blue")
lines(a3$it, a3$error, type = "l", xlab = "Iterations", ylab = "Relative Error",  col = "green")  
legend("topright", legend=c("Jacobi", "Parallel Jacobi", "Gauss-Seidel"),
       col=c("red", "blue", "green"), lty = 1, cex=0.8)
dev.off()


pdf(file="plot_time_2.pdf")
plot(a1$t, a1$error, type = "l", xlab = "Milliseconds", ylab = "Relative Error",  col = "red" )
lines(a2$t, a2$error, type = "l", xlab = "Milliseconds", ylab = "Relative Error",  col = "blue")
lines(a3$t, a3$error, type = "l", xlab = "Milliseconds", ylab = "Relative Error",  col = "green")  
legend("topright", legend=c("Jacobi", "Parallel Jacobi", "Gauss-Seidel"),
       col=c("red", "blue", "green"), lty = 1, cex=0.8)

dev.off()
