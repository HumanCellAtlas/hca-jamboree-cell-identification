
downsample=function(u,min_umis){
  all_op=1:nrow(u)
  base_tab=rep(0,nrow(u))
  names(base_tab)=all_op
  
  downsamp_one=function(x,n){
    tab=base_tab
    tab2=table(sample(rep(all_op,x),size = n,replace=F))
    tab[names(tab2)]=tab2
    return(tab)
  }
}