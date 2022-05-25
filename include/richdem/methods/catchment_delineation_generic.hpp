#ifndef _richdem_catchment_delineation_generic_hpp_
#define _richdem_catchment_delineation_generic_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/Array3D.hpp>
#include <richdem/common/logger.hpp>
#include <richdem/common/ProgressBar.hpp>
#include <richdem/methods/d8_methods.hpp> //only for sgn() function: TODO
#include <queue>

namespace richdem {

template<class T, class U>
std::queue<GridCell> queue_from_mask(
  const Array2D<U>  &mask,
  Array2D<T>              &upslope_cells
){
  std::queue<GridCell> expansion;
  // this cannot be parallelized, because it appends everythign to the same queue
  for (auto x=0;x<mask.width();x++){
    for (auto y=0;y<mask.height();y++){
      if(mask(x,y)!=mask.noData()){
        upslope_cells(x,y)=2;
        expansion.emplace(x,y);
      }
    }
  }

  return expansion;
}

template<class T>
std::queue<GridCell> queue_from_linepoints(
  int x0,
  int y0,
  int x1,
  int y1,
  Array2D<T>       &upslope_cells
){
  std::queue<GridCell> expansion;

  if(x0>x1){
    std::swap(x0,x1);
    std::swap(y0,y1);
  }

  //Modified Bresenham Line-Drawing Algorithm
  int deltax     = x1-x0;
  int deltay     = y1-y0;
  float error    = 0;
  float deltaerr = (float)deltay/(float)deltax;

  if (deltaerr<0)
    deltaerr = -deltaerr;

  RDLOG_MISC<<"Line slope is "<<deltaerr;
  int y=y0;
  for(int x=x0;x<=x1;x++){
    expansion.push(GridCell(x,y));
    upslope_cells(x,y)=2;
    error+=deltaerr;
    if (error>=0.5) {
      expansion.push(GridCell(x+1,y));
      upslope_cells(x+1,y) = 2;
      y                   += sgn(deltay);
      error               -= 1;
    }
  }
 
  return expansion;
}

//upslope_from_props
/**
 * @brief Calculates uplsope cells for a given set of input cells, given a
 * proportions matrix of flow. All cells that have flow into input cells will be added
 * 
 * @param[in]  &expansion      queue with starting cells
 * @param[in]  &elevation      DEM
 * @param[out] &upslope_cells  A grid of 1/2/NoData
 */
template<class T, class U>
void upslope_cells_props(
  std::queue<GridCell> &expansion,
  const Array3D<T>     &props,
  Array2D<U>           &upslope_cells
){
  //TODO: ALG NAME?
  RDLOG_PROGRESS<<"Setting up the upslope_cells matrix..."<<std::flush;
  // upslope_cells.resize(props.width(),props.height());
  // upslope_cells.setAll(0);
  // upslope_cells.setNoData(0);

  ProgressBar progress;

  progress.start(props.size());
  long int ccount=0;
  while(expansion.size()>0){
    GridCell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!props.inGrid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(upslope_cells.isNoData(c.x+dx[n],c.y+dy[n]) && props(c.x+dx[n],c.y+dy[n],d8_inverse[n])>0) {
        expansion.emplace(c.x+dx[n],c.y+dy[n]);
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop();
  RDLOG_MISC<<"Found "<<ccount<<" up-slope cells."; //TODO
}

//upslope_cells multiple flow implementation
/**
 * @brief Calculates uplsope cells for a given set of input cells, assuming
 * multiple flow. That is, all cells that are higher than the cell being
 * processed will be added, regardless of whether the current cell is the lowest
 * neighbour.
 * 
 * 
 * @param[in]  &expansion      queue with starting cells
 * @param[in]  &elevation      DEM
 * @param[out] &upslope_cells  A grid of 1/2/NoData
 */
template<class T, class U>
void upslope_cells_mflow(
  std::queue<GridCell> &expansion,
  const Array2D<T>     &elevation,
  Array2D<U>           &upslope_cells
){
  //TODO: ALG NAME?
  RDLOG_PROGRESS<<"Setting up the upslope_cells matrix..."<<std::flush;
  // upslope_cells.resize(elevation);
  // upslope_cells.setAll(0);
  // upslope_cells.setNoData(0);

  ProgressBar progress;

  progress.start(elevation.size());
  long int ccount=0;
  while(expansion.size()>0){
    GridCell c=expansion.front();
    expansion.pop();

    progress.update(ccount++);

    for(int n=1;n<=8;n++)
      if(!elevation.inGrid(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(elevation.isNoData(c.x+dx[n],c.y+dy[n]))
        continue;
      else if(upslope_cells.isNoData(c.x+dx[n],c.y+dy[n]) && elevation(c.x,c.y)<elevation(c.x+dx[n],c.y+dy[n])) {
        expansion.emplace(c.x+dx[n],c.y+dy[n]);
        upslope_cells(c.x+dx[n],c.y+dy[n])=1;
      }
  }
  RDLOG_TIME_USE<<"Succeeded in "<<progress.stop();
  RDLOG_MISC<<"Found "<<ccount<<" up-slope cells."; //TODO
}

}
#endif