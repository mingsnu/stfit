RMSE = function(y, ypred){
  sqrt(sum((y-ypred)^2)/length(y))
}
NMSE = function(y, ypred){
  sum((y-ypred)^2)/sum(y^2)
}
ARE = function(y, ypred){
  mean(abs(y - ypred)/y)
}