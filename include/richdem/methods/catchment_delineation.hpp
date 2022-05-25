#ifndef _richdem_catchment_delineation_hpp_
#define _richdem_catchment_delineation_hpp_

#include <richdem/common/Array2D.hpp>
#include <richdem/common/Array3D.hpp>
#include <richdem/methods/catchment_delineation_generic.hpp>

namespace richdem {

// TODO: do we have to resize  and initialize the results array?
template<class mask_t,               class result_t> void DC_mask_props(const Array2D<mask_t> &mask,    const Array3D<float> &props, Array2D<result_t> &upslope){std::queue<GridCell> expansion = queue_from_mask(mask,upslope)             ; upslope_cells_props(expansion, props, upslope);}
template<class mask_t, class elev_t, class result_t> void DC_mask_mflow(const Array2D<mask_t> &mask,    const Array2D<elev_t> &elev, Array2D<result_t> &upslope){std::queue<GridCell> expansion = queue_from_mask(mask,upslope)             ; upslope_cells_mflow(expansion, elev,  upslope);}
template<                            class result_t> void DC_line_props(int x0, int y0, int x1, int y1, const Array3D<float> &props, Array2D<result_t> &upslope){std::queue<GridCell> expansion = queue_from_linepoints(x0,y0,x1,y1,upslope); upslope_cells_props(expansion, props, upslope);}
template<              class elev_t, class result_t> void DC_line_mflow(int x0, int y0, int x1, int y1, const Array2D<elev_t> &elev, Array2D<result_t> &upslope){std::queue<GridCell> expansion = queue_from_linepoints(x0,y0,x1,y1,upslope); upslope_cells_mflow(expansion, elev,  upslope);}

}
#endif
