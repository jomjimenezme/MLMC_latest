#include "main.h"

void get_sigma_4D_8x8x8x8(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;

    g.num_colors[0] = 10;
  }
}

void get_sigma_3D_8x8x8(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;

    g.num_colors[0] = 8;
  }
}

void get_sigma_4D_64x32x32x32(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;

    g.num_colors[0] = 10;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 1;
    g.sigma[1] = 5;
    g.sigma[2] = 55;
    g.sigma[3] = 61;

    g.num_colors[0] = 16;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 1;
    g.sigma[1] = 8;
    g.sigma[2] = 12;
    g.sigma[3] = 18;

    g.num_colors[0] = 64;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 38;
    g.sigma[1] = 1;
    g.sigma[2] = 12;
    g.sigma[3] = 16;

    g.num_colors[0] = 128;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 3;
    g.sigma[1] = 20;
    g.sigma[2] = 48;
    g.sigma[3] = 50;

    g.num_colors[0] = 320;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 40;
    g.sigma[1] = 32;
    g.sigma[2] = 33;
    g.sigma[3] = 61;

    g.num_colors[0] = 512;
  }

}

void get_sigma_3D_32x32x32(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 2;
    g.sigma[3] = 3;

    g.num_colors[0] = 8;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 3;
    g.sigma[3] = 5;

    g.num_colors[0] = 16;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 6;
    g.sigma[3] = 9;

    g.num_colors[0] = 32;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 11;
    g.sigma[3] = 27;

    g.num_colors[0] = 88;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 8;
    g.sigma[3] = 44;

    g.num_colors[0] = 128;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 9;
    g.sigma[3] = 33;

    g.num_colors[0] = 176;
  }

  if(g.coloring_distance == 8){
    g.sigma[0] = 0;
    g.sigma[1] = 7;
    g.sigma[2] = 48;
    g.sigma[3] = 51;

    g.num_colors[0] = 272;
  }

  if(g.coloring_distance == 9){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 33;
    g.sigma[3] = 45;

    g.num_colors[0] = 352;
  }

}

void get_sigma_4D_128x64x64x64(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;

    g.num_colors[0] = 16;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 1;
    g.sigma[1] = 3;
    g.sigma[2] = 5;
    g.sigma[3] = 7;

    g.num_colors[0] = 16;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 1;
    g.sigma[1] = 4;
    g.sigma[2] = 10;
    g.sigma[3] = 17;

    g.num_colors[0] = 64;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 1;
    g.sigma[1] = 6;
    g.sigma[2] = 16;
    g.sigma[3] = 52;

    g.num_colors[0] = 128;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 6;
    g.sigma[1] = 1;
    g.sigma[2] = 48;
    g.sigma[3] = 120;

    g.num_colors[0] = 256;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 1;
    g.sigma[1] = 7;
    g.sigma[2] = 33;
    g.sigma[3] = 89;

    g.num_colors[0] = 448;
  }

  if(g.coloring_distance == 8){
    g.sigma[0] = 4;
    g.sigma[1] = 26;
    g.sigma[2] = 114;
    g.sigma[3] = 117;

    g.num_colors[0] = 832;
  }

}

void get_sigma_3D_64x64x64(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 3;
    g.sigma[3] = 5;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 2;
    g.sigma[3] = 3;

    g.num_colors[0] = 8;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 3;
    g.sigma[3] = 5;

    g.num_colors[0] = 16;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 4;
    g.sigma[3] = 10;

    g.num_colors[0] = 32;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 5;
    g.sigma[3] = 13;

    g.num_colors[0] = 64;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 6;
    g.sigma[3] = 26;

    g.num_colors[0] = 128;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 7;
    g.sigma[3] = 57;

    g.num_colors[0] = 152;
  }

  if(g.coloring_distance == 8){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 36;
    g.sigma[3] = 52;

    g.num_colors[0] = 256;
  }

  if(g.coloring_distance == 9){
    g.sigma[0] = 0;
    g.sigma[1] = 4;
    g.sigma[2] = 30;
    g.sigma[3] = 35;

    g.num_colors[0] = 320;
  }

}


void get_sigma_4D_192x96x96x96(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 1;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 1;
    g.sigma[1] = 2;
    g.sigma[2] = 3;
    g.sigma[3] = 4;

    g.num_colors[0] = 12;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 1;
    g.sigma[1] = 3;
    g.sigma[2] = 5;
    g.sigma[3] = 7;

    g.num_colors[0] = 16;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 1;
    g.sigma[1] = 4;
    g.sigma[2] = 11;
    g.sigma[3] = 17;

    g.num_colors[0] = 48;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 1;
    g.sigma[1] = 5;
    g.sigma[2] = 13;
    g.sigma[3] = 29;

    g.num_colors[0] = 96;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 1;
    g.sigma[1] = 6;
    g.sigma[2] = 50;
    g.sigma[3] = 76;

    g.num_colors[0] = 192;
  }

}

void get_sigma_3D_96x96x96(){

  if(g.coloring_distance == 1){
    g.sigma[0] = 0;
    g.sigma[1] = 1;
    g.sigma[2] = 1;
    g.sigma[3] = 1;

    g.num_colors[0] = 2;
  }

  if(g.coloring_distance == 2){
    g.sigma[0] = 0;
    g.sigma[1] = 5;
    g.sigma[2] = 1;
    g.sigma[3] = 6;

    g.num_colors[0] = 8;
  }

  if(g.coloring_distance == 3){
    g.sigma[0] = 0;
    g.sigma[1] = 7;
    g.sigma[2] = 3;
    g.sigma[3] = 1;

    g.num_colors[0] = 12;
  }

  if(g.coloring_distance == 4){
    g.sigma[0] = 0;
    g.sigma[1] = 28;
    g.sigma[2] = 31;
    g.sigma[3] = 19;

    g.num_colors[0] = 32;
  }

    if(g.coloring_distance == 5){
    g.sigma[0] = 0;
    g.sigma[1] = 41;
    g.sigma[2] = 23;
    g.sigma[3] = 13;

    g.num_colors[0] = 48;
  }

  if(g.coloring_distance == 6){
    g.sigma[0] = 0;
    g.sigma[1] = 21;
    g.sigma[2] = 65;
    g.sigma[3] = 30;

    g.num_colors[0] = 96;
  }

  if(g.coloring_distance == 7){
    g.sigma[0] = 0;
    g.sigma[1] = 42;
    g.sigma[2] = 49;
    g.sigma[3] = 117;

    g.num_colors[0] = 144;
  }

}
