#pragma once

#include <TCanvas.h>
#include <TH1.h>

static int HISTO_INDEX = 0; //current canvas number
static int CANVAS_INDEX = 0; //current canvas number

TCanvas * get_new_tailed_canvas(std::string title)
{
  auto c = new TCanvas;
  const int Nx = 4; //numbe of canvases on x size;
  const int Ny = 3; //number of canvases on y size;
  const int window_title_bar_ysize = 45;
  const int window_border_xsize = 5;
  const int panel_ysize = 40;
  const int display_size_x = 1920*2;
  const int display_size_y = 1080*2;
  const int canvas_width_x = display_size_x/Nx;
  const int canvas_width_y = (display_size_y-panel_ysize)/Ny;
  int canvas_pos_x = (CANVAS_INDEX % Nx)* canvas_width_x;
  int canvas_pos_y = ((CANVAS_INDEX/Nx)%Ny)* canvas_width_y;
  c->SetWindowPosition(canvas_pos_x,canvas_pos_y);
  c->SetWindowSize(canvas_width_x-window_border_xsize,canvas_width_y-window_title_bar_ysize);
  c->SetTitle(title.c_str());
  c->cd();
  CANVAS_INDEX++;
  return c;
}
