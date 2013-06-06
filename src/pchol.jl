function pchol(A,tile_sz)

# 1. Set up Data Structure
tile_sz2 =	tile_sz^2;
mat_sz   =	size(A,1);
bd_sz    =	tile_sz;
num_tile =	int(ceil(mat_sz/tile_sz));
upper = repmat(linspace(1,tile_sz,tile_sz)',tile_sz,1).>repmat(linspace(1,tile_sz,tile_sz),1,tile_sz);
tile_upper =	find(upper);
bd_upper   =	tile_upper;

@sync begin
#A_para is a cell, where each cell is the corresponding tile of matrix A
A_para = cell(num_tile,num_tile);
for i=1:num_tile-1
    for j=1:i
        A_para[i,j]=@spawn reshape(A[(i-1)*tile_sz+(1:tile_sz),(j-1)*tile_sz+(1:tile_sz)],1,tile_sz2);
    end
end

if tile_sz*num_tile>mat_sz
    bd_sz = mat_sz-(num_tile-1)*tile_sz;
    bd_upper=find(upper[1:bd_sz,1:bd_sz]);
    for j=1:num_tile-1
        A_para[num_tile,j]=@spawn reshape(A[(mat_sz-bd_sz+1):mat_sz,(j-1)*tile_sz+(1:tile_sz)],1,tile_sz*bd_sz);
    end
    A_para[num_tile,num_tile]=@spawn reshape(A[(mat_sz-bd_sz+1):mat_sz,(mat_sz-bd_sz+1):mat_sz],1,bd_sz*bd_sz);
else
    for j=1:num_tile
        A_para[num_tile,j]=@spawn reshape(A[(num_tile-1)*tile_sz+(1:tile_sz),(j-1)*tile_sz+(1:tile_sz)],1,tile_sz2);
    end
end
end

c_uplo   =	'L'; 
c_side   =	'R';
c_transt =	'T';
c_transn =	'N';
c_diag   =	'N';
c_alphap =	1.0;
c_alpham =	-1.0;
c_major  =	102;
c_info   =	0;


# some maths...
# A= LL'where L is lower diagnoal
# let L =[L11,0;L21,L22], then A=[L11*L11',L11*L21';L21*L11',L21*L21'+L22*L22']
# k=1;  2.1 solve L11 (with A11); 2.2 solve L21 (with L11,A21); 2.3/2.4 solve L22 (second iter)

# 2. Computation
# 2.1 DPOTRF


for k=1:num_tile-1
    akk = A_para[k,k];
    A_para[k,k] = @spawn dpotrf(c_uplo,tile_sz,fetch(akk),c_info)
    #fetch(A_para[k,k])[tile_upper]=0;

#    println(A_para[k,k])
#    println(fetch(A_para[k,k]))
#    wait(@spawnat A_para[k,k].where fetch(A_para[k,k])[tile_upper]=0)

    akk = fetch(A_para[k,k]);
    akk[tile_upper] = 0;
    A_para[k,k] = akk; 
#    println(fetch(A_para[k,k]))
#println(k)
    # 2.2 DTRSM
    for m=k+1:num_tile-1
        akk = A_para[k,k]; 
	amk = A_para[m,k];
        A_para[m,k] = @spawn dtrsm(c_side,c_uplo,c_transt,c_diag,tile_sz,tile_sz,c_alphap,
		fetch(akk),fetch(amk))   
    end
    
    akk = A_para[k,k];
    antk = A_para[num_tile,k];
    A_para[num_tile,k] = @spawn dtrsm(c_side,c_uplo,c_transt,c_diag,bd_sz,tile_sz,c_alphap,
	fetch(akk),fetch(antk))

    # 2.3 DSYRK

    for n=k+1:num_tile-1
	ann = A_para[n,n];
	ank = A_para[n,k];
	A_para[n,n] = @spawn dsyrk(c_uplo,c_transn,tile_sz,tile_sz,c_alpham,fetch(ank),
		c_alphap,fetch(ann))

        # 2.4 DGEMM

        for m=n+1:num_tile-1
	    ank = A_para[n,k];
	    amk = A_para[m,k];
	    amn = A_para[m,n];
	    A_para[m,n] = @spawn dgemm(c_transn,c_transt,tile_sz,tile_sz,c_alpham,
		fetch(amk),fetch(ank),c_alphap,fetch(amn)) 
        end

	ank = A_para[n,k];
	antn = A_para[num_tile,n];
	antk = A_para[num_tile,k];
	A_para[num_tile,n] = @spawn dgemm(c_transn,c_transt,bd_sz,tile_sz,c_alpham,
		fetch(antk),fetch(ank),c_alphap,fetch(antn))

    end

    antk = A_para[num_tile,k];
    antnt = A_para[num_tile,num_tile];
    A_para[num_tile,num_tile] = @spawn dsyrk(c_uplo,c_transn,bd_sz,tile_sz,c_alpham,
	fetch(antk), c_alphap, fetch(antnt))

end

# last diagonal block
antnt = A_para[num_tile,num_tile];
A_para[num_tile,num_tile] = @spawn dpotrf(c_uplo,bd_sz,fetch(antnt),c_info)

#println(A_para[num_tile,num_tile])

#wait(@spawnat A_para[num_tile,num_tile].where fetch(A_para[num_tile,num_tile])[bd_upper]=0);

antnt = fetch(A_para[num_tile,num_tile]);
antnt[bd_upper] = 0;
A_para[num_tile,num_tile] = antnt;
#A_para[num_tile,num_tile][bd_upper]=0;

A_out=zeros(size(A))
for i=1:num_tile-1
    for j=1:i#num_tile
        A_out[(i-1)*tile_sz+(1:tile_sz),(j-1)*tile_sz+(1:tile_sz)]=reshape(fetch(A_para[i,j]),tile_sz,tile_sz);
    end
end
if bd_sz!=tile_sz
    for j=1:num_tile-1
        A_out[(mat_sz-bd_sz+1):mat_sz,(j-1)*tile_sz+(1:tile_sz)]=reshape(fetch(A_para[num_tile,j]),bd_sz,tile_sz);
    end
    A_out[(mat_sz-bd_sz+1):mat_sz,(mat_sz-bd_sz+1):mat_sz]=reshape(fetch(A_para[num_tile,num_tile]),bd_sz,bd_sz);
else
    for j=1:num_tile
        A_out[(num_tile-1)*tile_sz+(1:tile_sz),(j-1)*tile_sz+(1:tile_sz)]=reshape(fetch(A_para[num_tile,j]),tile_sz,tile_sz);
    end
end

A_out

end

#http://www.netlib.org/clapack/clapack-3.2.1-CMAKE/SRC/VARIANTS/cholesky/TOP/dpotrf.c
#http://www.netlib.org/clapack/cblas/dtrsm.c
#http://www.netlib.org/clapack/cblas/dsyrk.c
#http://www.netlib.org/clapack/cblas/dgemm.c
