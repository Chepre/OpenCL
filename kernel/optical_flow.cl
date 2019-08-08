#defineIx(im1,i,j,W,H) ((P(im1, (i)+1, (j), W, H) -P(im1, (i)-1, (j), W, H)) >> 1)
#defineIy(im1,i,j,W,H) ((P(im1, (i), (j)+1, W, H) -P(im1, (i), (j)-1, W, H)) >> 1)

kernel voidoptical_flow_for_images(constant uchar* restrict im1, global uchar*restrict im2, intW_in, intH_in, global uchar* restrict out, constant int* restrict colorwheel, intncols)
{
	while(j < H -WINDOW_SIZE) 
	{
		float	G_inv[4]; 
                float	mu[2]; 
		int	G[4]; 
		float	b_k[2];
		int	wj= -WINDOW_SIZE, wi= -WINDOW_SIZE;

		while(wj<= WINDOW_SIZE) 
		{
			// Image difference. Formula #30u
			char	im2_val = P(im2, (int)(i+wi), (int)(j+wj), W, H);
			int	deltaIk= P(im1, i+wi, j+wj, W, H) -im2_val;
			int	cIx= Ix(im1, i+wi,j+wj,W,H);
			int	cIy= Iy(im1, i+wi,j+wj,W,H);
			
			G[0] += cIx* cIx;  
			G[1] += cIx* cIy;
			G[2] += cIx* cIy;  
			G[3] += cIy* cIy;
			b_k[0] += deltaIk* cIx;   
			b_k[1] += deltaIk* cIy;

			wj= (wi>WINDOW_SIZE) ? (wj+1) : wj;
			wi= (wi>WINDOW_SIZE) ? -WINDOW_SIZE : (wi+1);
		}    

		grad_invertible= get_matrix_inv(G, G_inv);

		// Guess for next iteration. Formulas #28 & #31
		mu[0] = mu[0] + G_inv[0] * b_k[0] + G_inv[1] * b_k[1];
		mu[1] = mu[1] + G_inv[2] * b_k[0] + G_inv[3] * b_k[1];
		float	fx = mu[0] + 0; 

		//2 * P(prev_g, 2*(i/2),   j/2, W, H/2);
		float	fy = mu[1] + 0; 
		//2 * P(prev_g, 2*(i/2)+1, j/2, W, H/2);
		global uchar *ptr = &(P(out, 3*i, j, 3*W, H));

		if(grad_invertible) 
		{
			computeColor(colorwheel, ncols, fx, fy, ptr);  
		}

		j = (i+1)>=(W-WINDOW_SIZE) ? (j+1) : j;
		i = (i+1)>=(W-WINDOW_SIZE) ? WINDOW_SIZE : (i+1);
	} 
}
