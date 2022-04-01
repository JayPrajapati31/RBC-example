/*Target profiles with balde_len = 0.5 & t in [0,1]*/

xy[2*i+1] = 15*xy[2*i]*(xy[2*i]-0.125)*(xy[2*i]-0.25) + 0.4;
xy[2*i+1] = 150*xy[2*i]*(xy[2*i]-0.1*block_size)*(xy[2*i]-0.9*block_size)*(xy[2*i]-0.25) + 0.4;
xy[2*i+1] = 300*xy[2*i]*(xy[2*i]-0.4*block_size)*(xy[2*i]-0.6*block_size)*(xy[2*i]-0.25) + 0.4;
xy[2*i+1] = 0.1*xy[2*i]*sin(8*M_PI*xy[2*i]) + 0.3;
xy[2*i+1] = 0.1*xy[2*i]*sin(14*M_PI*xy[2*i]) + 0.4;
xy[2*i+1] = 0.01*sin(16*M_PI*xy[2*i]) + 0.5;