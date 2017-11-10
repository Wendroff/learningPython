Module OUTFLOW
    implicit none
    real*8,allocatable :: S(:),DENe(:),Ue(:),We(:),Te(:) !无量纲化的结果
    real*8             :: S0 !Fai(1,:)的实际无量纲位置，计算时从S上统一减掉，最后又加回来
                             !之所以这样，是因为自启动的初始条件必须x=0
    !logical            :: Compress_flag
    
    
    
end module