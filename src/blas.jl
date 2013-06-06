function dgeqrt(m,n,nb,A,lda,T,ldt,work,info)

	ccall(dlsym(libBLAS, :dgeqrt_),Void,
	(Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32},
	Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}),
	&m, &n, &nb, A, &lda, T, &ldt, work, &info)

	(A,T,work)

end

function dtpqrt(m,n,l,nb,A,lda,B,ldb,T,ldt,work,info)

	ccall(dlsym(libBLAS, :dtpqrt_),Void,
	(Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
	Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, 
	Ptr{Float64}, Ptr{Int32}),
	&m, &n, &l, &nb, A, &lda, B, &ldb, T, &ldt, work, &info)

	(A,B,T,work)

end


function dormqr(side,trans,m,n,k,A,lda,T,C,ldc,work,lwork,info)

	ccall(dlsym(libBLAS, :dormqr_),Void,
	(Ptr{Uint8}, Ptr{Uint8}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
	Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32},
	Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
	&side,&trans,&m,&n,&k,A,&lda,T,C,&ldc,work,&lwork,&info)

	(C,work)

end

function dtpmqrt(side,trans,m,n,k,l,nb,V,ldv,T,ldt,A,lda,B,ldb,work,info)

        ccall(dlsym(libBLAS, :dtpmqrt_),Void,
	(Ptr{Uint8}, Ptr{Uint8}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
	Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64},
	Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32},
	Ptr{Float64}, Ptr{Int32}),
	&side,&trans,&m,&n,&k,&l,&nb,V,&ldv,T,&ldt,A,&lda,B,&ldb,work,&info)

	(A,B,work)

end

function dpotrf(c_uplo,tile_sz,A,c_info)

	ccall(dlsym(libBLAS, :dpotrf_),Void,
	(Ptr{Uint8}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}),
	&c_uplo, &tile_sz, A, &tile_sz, &c_info)

	A

end

function dtrsm(c_side,c_uplo,c_transt,c_diag,bd_sz,tile_sz,c_alphap,A1,A2)

	ccall(dlsym(libBLAS, :dtrsm_),Void,
	(Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32},Ptr{Int32},
	Ptr{Float64},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Int32}),
	&c_side, &c_uplo,&c_transt,&c_diag,&bd_sz,&tile_sz,&c_alphap,
	A1,&tile_sz,A2,&bd_sz)

	A2

end

function dsyrk(c_uplo,c_transn,bd_sz,tile_sz,c_alpham,A1,c_alphap,A2)

	ccall(dlsym(libBLAS, :dsyrk_),Void,
	(Ptr{Uint8},Ptr{Uint8},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},
	Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Int32}),
	&c_uplo,&c_transn,&bd_sz,&tile_sz,&c_alpham,A1,
	&bd_sz,&c_alphap,A2, &bd_sz)

	A2

end


function dgemm(c_transn,c_transt,bd_sz,tile_sz,c_alpham,A1,A2,c_alphap,A3)

	ccall(dlsym(libBLAS, :dgemm_),Void,
	(Ptr{Uint8},Ptr{Uint8},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},
	Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Float64},
	Ptr{Float64},Ptr{Int32}),
	&c_transn,&c_transt,&bd_sz, &tile_sz,&tile_sz,&c_alpham,A1,&bd_sz,
	A2,&tile_sz,&c_alphap,A3,&bd_sz)

	A3

end
