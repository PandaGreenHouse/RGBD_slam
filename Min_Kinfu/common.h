#pragma once
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/common/common.h>
#include <pcl/pcl_config.h>
#include <pcl/common/transforms.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/morphological_filter.h>
#include <pcl/features/moment_of_inertia_estimation.h>
#include <pcl/common/pca.h>
#include <pcl/segmentation/supervoxel_clustering.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/visualization/point_cloud_color_handlers.h>
#include <pcl/io/ply_io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/kdtree/kdtree.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/extract_clusters.h>
//#include "Simplify.h"
void refineMesh(pcl::PolygonMesh &mesh);
void clipping(pcl::PolygonMesh &mesh_in, float fRadius);
void removeFloor(pcl::PolygonMesh &mesh_in, float plane_tolerance);
void extractObject(pcl::PolygonMesh &mesh_in, float clip_radius, float plane_tolerance);
void decimateMesh(pcl::PolygonMesh& mesh_in, pcl::PolygonMesh& mesh_out);
void drawScene(Eigen::Affine3f &camPose, int x, int y, int width, int height, std::string &caption, boost::shared_ptr<pcl::PolygonMesh> mesh_ptr);
