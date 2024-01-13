# code for regression analysis

# artificial data
x = c(100, 110, 120, 130, 140, 150, 160, 170, 180, 190)
y = c( 45,  52,  54,  63,  62,  68,  75,  76,  92,  88)
plot(x, y)

# fiber data
x = c(1, 2, 3, 4, 5, 6, 7)
y = c(3, 2, 5, 6, 4, 8, 9)
fit = lm(y ~ x)
print(fit)

# moisture data
x = c(46, 53, 29, 61, 36, 39, 47, 49, 52, 38, 55, 32, 57, 54, 44)
y = c(12, 15,  7, 17, 10, 11, 11, 12, 14,  9, 16,  8, 18, 14, 12)
moisture = lm(y ~ x)
print(moisture)
plot(x, y)
abline(moisture)

# consumption data
x = c(  45, 50,   55, 60,   65,   70,   75)
y = c(24.2, 25, 23.3, 22, 21.5, 20.6, 19.8)
miles = lm(y ~ x)
print(miles)
print(summary(miles))

# father-son data
x = c(  60,   62, 64,   65,   66,   67,   68,   70,   72, 74)
y = c(63.6, 65.2, 66, 65.5, 66.9, 67.1, 67.4, 68.3, 70.1, 70)
father_son = lm(y ~ x)
print("Father-Son Data")
print(father_son)
print(summary(father_son))

