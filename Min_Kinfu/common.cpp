#include "common.h"
#include "Simplify.h"
using namespace pcl;
bool loop = true;
pcl::visualization::PCLVisualizer::Ptr _viewer;
boost::shared_ptr<pcl::PolygonMesh> _mesh_ptr;
std::vector<std::vector<int>> _adjancencies;
pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_total(new pcl::PointCloud<pcl::PointXYZRGB>);
std::vector<pcl::PointCloud<pcl::PointXYZRGB>::Ptr> cloud_clusters;
std::vector<std::pair<int, int>> clusters;
int cluster_no = 0;
void refineMesh(PolygonMesh &mesh)
{
	PointCloud<PointXYZ> cloud, cloud_refined;
	pcl::fromPCLPointCloud2(mesh.cloud, cloud);
	std::vector<int> labels;
	labels.resize(cloud.points.size(), -1);
	size_t count = 0;
	for (size_t i = 0; i < labels.size(); ++i)
	{
		if (labels[i] > -1)
			continue;
		labels[i] = count;
		cloud_refined.points.push_back(cloud.points[i]);
		Eigen::Vector3f v1 = cloud.points[i].getVector3fMap();
		for (size_t j = 0; j < cloud.points.size(); ++j)
		{
			if (labels[j] == -1)
			{
				Eigen::Vector3f v2 = cloud.points[j].getVector3fMap();
				if (v1 == v2)
					labels[j] = count;
			}
		}
		++count;
	}
	for (size_t i = 0; i < mesh.polygons.size(); ++i)
	{
		Vertices vts = mesh.polygons.at(i);
		vts.vertices[0] = labels[vts.vertices[0]];
		vts.vertices[1] = labels[vts.vertices[1]];
		vts.vertices[2] = labels[vts.vertices[2]];
	}
	pcl::toPCLPointCloud2(cloud_refined, mesh.cloud);
}

void generatingAdjancencies()
{
	_adjancencies.resize(_mesh_ptr->polygons.size());
	for (size_t i = 0; i < _mesh_ptr->polygons.size(); ++i)
	{
		Vertices verts1 = _mesh_ptr->polygons.at(i);
		_adjancencies[i].resize(0);
		for (size_t j = 0; j < _mesh_ptr->polygons.size(); ++j)
		{
			if (i == j)
				continue;
			Vertices verts2 = _mesh_ptr->polygons.at(j);
			int count = 0;
			for (int k = 0; k < 3; ++k)
			{
				for (int n = 0; n < 3; ++n)
				{
					if (verts1.vertices[k] == verts2.vertices[n])
					{
						++count;
						break;
					}
				}
				if (count == 2)
				{
					_adjancencies[i].push_back(j);
					break;
				}
			}
		}
	}
}
/*Nov, 30, 2022*/
void clipping(PointCloud<pcl::PointXYZ> &cloud_in, float clip_radius)
{
	pcl::PointCloud<pcl::PointXYZ> cloud_clipped;
	Eigen::Vector3f vCenter(0.0f, 0.0f, 0.0f);
	for (size_t i = 0; i < cloud_in.points.size(); ++i)
	{
		vCenter = (i*vCenter + cloud_in.points[i].getVector3fMap()) / (i + 1);
	}
	cloud_clipped.points.resize(0);
	for (size_t i = 0; i < cloud_in.points.size(); ++i)
	{
		Eigen::Vector3f vec = cloud_in.points[i].getVector3fMap() - vCenter;
		if (vec.norm() < clip_radius)
		{
			cloud_clipped.points.push_back(cloud_in.points[i]);
		}
	}
	cloud_in.points.swap(cloud_clipped.points);
}

void clipping(pcl::PolygonMesh &mesh_in, float radius)
{
	Eigen::Vector3f vCenter(0.0f, 0.0f, 0.0f);
	PointCloud<PointXYZ> cloud_in, cloud_clipped;
	std::vector<int> indices;
	pcl::fromPCLPointCloud2(mesh_in.cloud, cloud_in);
	indices.resize(0);
	cloud_clipped.points.resize(0);
	/*for (int i = 0; i < cloud_in.points.size(); ++i)
	{
		vCenter = (i*vCenter + cloud_in.points[i].getVector3fMap()) / (i + 1);
	}*/
	vCenter[0] = 0.5f*radius;
	vCenter[1] = 0.5f*radius;
	vCenter[2] = 0.5f*radius;
	std::vector<int> labels;
	labels.resize(cloud_in.points.size());
	size_t count = 0;
	for (size_t i = 0; i < cloud_in.points.size(); ++i)
	{
		labels[i] = -1;
		Eigen::Vector3f vec = cloud_in.points[i].getVector3fMap() - vCenter;
		if (vec.norm() > radius)
			continue;
		indices.push_back(i);
		cloud_clipped.points.push_back(cloud_in.points[i]);
		labels[i] = count;
		++count;
	}
	std::vector<Vertices> polygons_clipped;
	polygons_clipped.resize(0);
	for (size_t i = 0; i < mesh_in.polygons.size(); ++i)
	{
		Vertices t;
		t.vertices.resize(3);
		bool flg = true;
		for (int j = 0; j < 3; ++j)
		{
			int id = mesh_in.polygons.at(i).vertices[j];
			if (labels[id] == -1)
			{
				flg = false;
				break;
			}
			t.vertices[j] = labels[id];
		}
		if (flg == false)
			continue;
		polygons_clipped.push_back(t);
	}
	pcl::toPCLPointCloud2(cloud_clipped, mesh_in.cloud);
	mesh_in.polygons.swap(polygons_clipped);
}

void removeFloor(pcl::PolygonMesh &mesh_in, float plane_tolerance)
{
	PointCloud<PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>());
	pcl::fromPCLPointCloud2(mesh_in.cloud, *cloud_in);
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients());
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices());
	// Create the segmentation object
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	// Optional
	seg.setOptimizeCoefficients(true);
	// Mandatory
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setMaxIterations(200);
	//float plane_tolerance = 0.05f;
	seg.setDistanceThreshold(plane_tolerance);//(0.05);
	seg.setInputCloud(cloud_in);
	seg.segment(*inliers, *coefficients);
	std::vector<int> map;
	map.resize(cloud_in->points.size(),1);
	for (size_t i = 0; i < inliers->indices.size(); ++i)
	{
		map[inliers->indices[i]] = -1;
	}
	size_t count = 0;
	PointCloud<PointXYZ> cloud_final;
	cloud_final.points.resize(0);
	for (size_t i = 0; i < map.size(); ++i)
	{
		if (map[i] == -1)
			continue;
		map[i] = count;
		cloud_final.push_back(cloud_in->points[i]);
		++count;
	}
	std::vector<Vertices> polygons;
	polygons.resize(0);
	for (size_t i = 0; i < mesh_in.polygons.size(); ++i)
	{
		Vertices vts = mesh_in.polygons.at(i);
		bool flg = true;
		for (int j = 0; j < 3; ++j)
		{
			if (map[vts.vertices[j]] == -1)
			{
				flg = false;
				break;
			}
		}
		if (flg == false)
			continue;
		Vertices t;
		t.vertices.resize(3);
		for (int j = 0; j < 3; ++j)
		{
			t.vertices[j] = map[vts.vertices[j]];
		}
		polygons.push_back(t);
	}
	PointCloud<PointXYZRGB> cloud_rgb;
	for (size_t i = 0; i < cloud_final.points.size(); ++i)
	{
		PointXYZRGB p;
		p.x = cloud_final.points[i].x;
		p.y = cloud_final.points[i].y;
		p.z = cloud_final.points[i].z;
		p.r = 150;
		p.g = 150;
		p.b = 150;
		cloud_rgb.points.push_back(p);
	}
	pcl::toPCLPointCloud2(cloud_rgb, mesh_in.cloud);
	mesh_in.polygons.swap(polygons);
}

void extractObject(PolygonMesh &mesh_in, float clip_radius, float plane_tolerance)
{
	//clipping(mesh_in, clip_radius);
	removeFloor(mesh_in, plane_tolerance);
}

void decimateMesh(pcl::PolygonMesh& mesh_in, pcl::PolygonMesh& mesh_out)
{
	pcl::PointCloud<pcl::PointXYZRGB> cloud_in;
	pcl::fromPCLPointCloud2(mesh_in.cloud, cloud_in);
	std::vector<Eigen::Vector3f> vertices_in;
	std::vector<int> indices_in;
	vertices_in.resize(0);
	indices_in.resize(0);
	for (int i = 0; i < cloud_in.points.size(); ++i)
	{
		Eigen::Vector3f v = cloud_in.points.at(i).getVector3fMap();
		vertices_in.push_back(v);
	}
	for (int i = 0; i < mesh_in.polygons.size(); ++i)
	{
		pcl::Vertices vertices = mesh_in.polygons[i];
		int id = vertices.vertices[0];
		indices_in.push_back(id);
		id = vertices.vertices[1];
		indices_in.push_back(id);
		id = vertices.vertices[2];
		indices_in.push_back(id);
	}
	Simplify::importMesh(vertices_in, indices_in);
	int target_count = mesh_in.polygons.size() / 10;//10,5;
	double agressiveness = 7;
	Simplify::simplify_mesh(target_count, agressiveness, true);
	std::vector<Eigen::Vector3f> vertices_decimated;
	std::vector<int> indices_decimated;
	Simplify::exportMesh(vertices_decimated, indices_decimated);
	pcl::PointCloud<pcl::PointXYZ> cloud_decimated;
	for (int i = 0; i < vertices_decimated.size(); ++i)
	{
		pcl::PointXYZ p;
		p.x = vertices_decimated[i][0];
		p.y = vertices_decimated[i][1];
		p.z = vertices_decimated[i][2];
		cloud_decimated.push_back(p);
	}
	pcl::toPCLPointCloud2(cloud_decimated, mesh_out.cloud);
	for (int i = 0; i < indices_decimated.size() / 3; ++i)
	{
		pcl::Vertices vertices;
		vertices.vertices.resize(3);
		vertices.vertices[0] = indices_decimated[3 * i];
		vertices.vertices[1] = indices_decimated[3 * i + 1];
		vertices.vertices[2] = indices_decimated[3 * i + 2];
		mesh_out.polygons.push_back(vertices);
	}
}

bool intersectTri(Eigen::Vector3f &v0, Eigen::Vector3f &v1, Eigen::Vector3f &v2,
	Eigen::Vector3f &vRayPos, Eigen::Vector3f &vRayDir, float &u, float &v, float &dist, Eigen::Vector3f &interPos)
{
	Eigen::Vector3f e1 = v1 - v0;
	Eigen::Vector3f e2 = v2 - v0;
	Eigen::Vector3f q = vRayDir.cross(e2);
	float a = e1.dot(q);//D3DXVec3Dot(&e1,&q);//e1.dot(q);
	Eigen::Vector3f s = vRayPos - v0;
	Eigen::Vector3f r = s.cross(e1);
	//D3DXVec3Cross(&r,&s,&e1);//s.cross(e1);
	// Barycentric vertex weights
	u = s.dot(q) / a;
	v = vRayDir.dot(r) / a;
	float w = 1.0f - (u + v);//weight[0] = 1.0f - (weight[1] + weight[2]);
	dist = e2.dot(r) / a;
	static const float epsilon = 1e-7f;
	static const float epsilon2 = 1e-10f;
	if ((a <= epsilon) || (u < -epsilon2) ||
		(v < -epsilon2) || (w < -epsilon2) ||
		(dist <= 0.0f)) {
		// The ray is nearly parallel to the triangle, or the
		// intersection lies outside the triangle or behind
		// the ray origin: "infinite" distance until intersection.
		return false;
	}
	else {
		interPos = v0 + u*e1 + v*e2;
		return true;
	}
}

int getIntersection(Eigen::Vector3f &vRayPos, Eigen::Vector3f &vRayDir, pcl::PolygonMesh &mesh, Eigen::Vector3f &intersection)
{
	bool bIntersect = false;
	pcl::PointCloud<pcl::PointXYZRGB> cloud_in;
	pcl::fromPCLPointCloud2(mesh.cloud, cloud_in);
	float u, v;
	Eigen::Vector3f pts[3];
	float minDist = std::numeric_limits<float>::max();
	int pickId = -1;
	for (size_t i = 0; i < mesh.polygons.size(); ++i)
	{
		//int indices[3];
		Eigen::Vector3f points[3];
		for (int j = 0; j < 3; ++j)
		{
			int id = mesh.polygons[i].vertices[j];
			points[j] = cloud_in.points[id].getVector3fMap();
		}
		float u1, v1, dist;
		Eigen::Vector3f interPos;
		if (intersectTri(points[0], points[1], points[2], vRayPos, vRayDir, u1, v1, dist, interPos))
		{
			if (minDist > dist)
			{
				for (int k = 0; k < 3; ++k)
				{
					pts[k] = points[k];
				}
				u = u1;
				v = v1;
				pickId = i;
				minDist = dist;
				intersection = interPos;
				bIntersect = true;
			}
		}
	}
	return pickId;
}

int pickingMesh(int px, int py, pcl::visualization::PCLVisualizer::Ptr viewer, pcl::PolygonMesh &mesh)
{
	std::vector<pcl::visualization::Camera> cams;
	viewer->getCameras(cams);
	int width = cams[0].window_size[0];
	int height = cams[0].window_size[1];
	float cx = (float)width / 2;
	float cy = (float)height / 2;
	float y_max = tan(cams[0].fovy / 2);
	float x_max = y_max*(float)width / height;
	float x = (float)px - cx;
	float y = (float)py - cy;
	Eigen::Vector3f p;
	float w = 0.5f*(float)width;
	float h = 0.5f*(float)height;
	p[0] = -x_max*(float)x / w;
	p[1] = y_max*(float)y / h;
	p[2] = 1.0f;
	Eigen::Affine3f pose = viewer->getViewerPose();
	Eigen::Vector3f vRayDir = pose.linear()*p;
	Eigen::Vector3f vRayPos = pose.translation();
	vRayDir.normalize();
	Eigen::Vector3f intersection;
	int faceId = getIntersection(vRayPos, vRayDir, mesh, intersection);
	if ( faceId > -1)
	{
		PointCloud<PointXYZRGB> cloud;
		pcl::fromPCLPointCloud2(mesh.cloud, cloud);
		Vertices vts = mesh.polygons.at(faceId);
		for (int i = 0; i < 3; ++i)
		{
			cloud.points[vts.vertices[i]].r = 255;
			cloud.points[vts.vertices[i]].g = 0;
			cloud.points[vts.vertices[i]].b = 0;
		}
		pcl::toPCLPointCloud2(cloud, mesh.cloud);
	}
	return faceId;
}

void segmentMesh(Eigen::Vector3f &vRayPos, Eigen::Vector3f &vRayDir, PolygonMesh &mesh)
{
	Eigen::Vector3f intersection;
	int facePicked = getIntersection(vRayPos, vRayDir, mesh, intersection);
	if ( facePicked > -1)
	{
		std::vector<int> labels;
		labels.resize(mesh.polygons.size(), -1);
		int count = 0;
		std::vector<int> seg;
		seg.push_back(facePicked);
		while (count < seg.size())
		{
			int t = seg[count];
			for (int i = 0; i < _adjancencies[t].size(); ++i)
			{
				if (labels[_adjancencies[t][i]] == -1)
				{
					labels[_adjancencies[t][i]] = 1;
					seg.push_back(_adjancencies[t][i]);
				}
			}
			++count;
		}
		PointCloud<PointXYZRGB> cloud;
		pcl::fromPCLPointCloud2(mesh.cloud, cloud);
		for (size_t i = 0; i < seg.size(); ++i)
		{
			Vertices vts = mesh.polygons[seg[i]];
			for (int j = 0; j < 3; ++j)
			{
				cloud.points[vts.vertices[j]].r = 0;
				cloud.points[vts.vertices[j]].g = 255;
				cloud.points[vts.vertices[j]].b = 0;
			}
		}
		pcl::toPCLPointCloud2(cloud, mesh.cloud);
	}
}

void segmentMesh(int facePicked, PolygonMesh &mesh)
{
	std::vector<int> labels;
	labels.resize(mesh.polygons.size(), -1);
	int count = 0;
	std::vector<int> seg;
	seg.push_back(facePicked);
	while (count < seg.size())
	{
		int t = seg[count];
		for (int i = 0; i < _adjancencies[t].size(); ++i)
		{
			if (labels[_adjancencies[t][i]] == -1)
			{
				labels[_adjancencies[t][i]] = 1;
				seg.push_back(_adjancencies[t][i]);
			}
		}
		++count;
	}
	PointCloud<PointXYZRGB> cloud;
	pcl::fromPCLPointCloud2(mesh.cloud, cloud);
	for (size_t i = 0; i < seg.size(); ++i)
	{
		Vertices vts = mesh.polygons[seg[i]];
		for (int j = 0; j < 3; ++j)
		{
			cloud.points[vts.vertices[j]].r = 0;
			cloud.points[vts.vertices[j]].g = 255;
			cloud.points[vts.vertices[j]].b = 0;
		}
	}
	pcl::toPCLPointCloud2(cloud, mesh.cloud);
}

void process_key_down(const pcl::visualization::KeyboardEvent& event, void* v)
{
	if ((event.getKeyCode() == 'a' || event.getKeyCode() == 'A') && event.keyUp() == true)
	{
		Eigen::Affine3f pose = _viewer->getViewerPose();
		Eigen::Vector3f vRayDir(0.0f, 1.0f, 0.0f);
		vRayDir = pose.linear()*vRayDir;
		Eigen::Vector3f pos = pose.translation();
		segmentMesh(pos, vRayDir, *_mesh_ptr);
	}
	if ((event.getKeyCode() == 's' || event.getKeyCode() == 'S') && event.keyUp() == true)
	{
		++cluster_no;
		if (cluster_no > clusters.size())
			cluster_no = 0;
	}
	if ((event.getKeyCode() == 'q' || event.getKeyCode() == 'Q'))
	{
		loop = false;
	}
	if ((event.getKeyCode() == 'c' || event.getKeyCode() == 'C') && event.keyUp() == true)
	{
		int no = cluster_no;
		++cluster_no;
		if (cluster_no > cloud_clusters.size())
			cluster_no = 0;
		cloud_clusters[no]->empty();
		cloud_clusters.erase(cloud_clusters.begin() + no);
	}
}

void process_mouse(const pcl::visualization::MouseEvent& event, void* v)
{
	if (event.getButton() == pcl::visualization::MouseEvent::RightButton)
	{
		if (event.getType() == pcl::visualization::MouseEvent::MouseButtonPress)
		{
			int px = event.getX();
			int py = event.getY();
			int facePicked = pickingMesh(px, py, _viewer, *_mesh_ptr);
			segmentMesh(facePicked, *_mesh_ptr);
		}
	}
}

void setViewerPose(pcl::visualization::PCLVisualizer::Ptr viewer, Eigen::Affine3f &pose)
{
	Eigen::Vector3f pos_vector = pose.translation();
	Eigen::Vector3f vLookAt = pos_vector + pose.linear()*Eigen::Vector3f::UnitZ();
	Eigen::Vector3f up_vec = pose.linear()*(-Eigen::Vector3f::UnitY());
	viewer->setCameraPosition(pos_vector[0], pos_vector[1], pos_vector[2], vLookAt[0], vLookAt[1], vLookAt[2], up_vec[0], up_vec[1], up_vec[2]);
}

pcl::visualization::PCLVisualizer::Ptr createViewer(Eigen::Affine3f &camPose, int x, int y, int width, int height, std::string &caption)
{
	pcl::visualization::PCLVisualizer::Ptr viewer = pcl::visualization::PCLVisualizer::Ptr(new pcl::visualization::PCLVisualizer(caption));
	viewer->setWindowName(caption);
	viewer->setBackgroundColor(0, 0, 0.15);
	viewer->addCoordinateSystem(1.0, "world");
	viewer->initCameraParameters();
	viewer->setPosition(x, y);
	viewer->setSize(width, height);
	viewer->setCameraClipDistances(0.0, 10.0);
	viewer->registerKeyboardCallback(&process_key_down);
	viewer->registerMouseCallback(&process_mouse);
	setViewerPose(viewer, camPose);
	return viewer;
}

void segmentPointCloud(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in, PointCloud<PointXYZRGB> &cloud_out, std::vector<std::pair<int, int>> &clusters)
{
	std::vector<pcl::RGB> colors;
	colors.resize(8);
	colors[0].r = 255; colors[0].g = 0; colors[0].b = 0;
	colors[1].r = 0; colors[1].g = 255; colors[1].b = 0;
	colors[2].r = 0; colors[2].g = 0; colors[2].b = 255;
	colors[3].r = 255; colors[3].g = 255; colors[3].b = 0;
	colors[4].r = 255; colors[4].g = 0; colors[4].b = 255;
	colors[5].r = 0; colors[5].g = 255; colors[5].b = 255;
	colors[6].r = 200; colors[6].g = 150; colors[6].b = 0;
	colors[7].r = 200; colors[7].g = 150; colors[7].b = 200;
	// Create the filtering object: downsample the dataset using a leaf size of 1cm
	pcl::VoxelGrid<pcl::PointXYZ> vg;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>), cloud_f(new pcl::PointCloud<pcl::PointXYZ>);
	vg.setInputCloud(cloud_in);
	vg.setLeafSize(0.01f, 0.01f, 0.01f);
	vg.filter(*cloud_filtered);
	// Create the segmentation object for the planar model and set all the parameters
	pcl::SACSegmentation<pcl::PointXYZ> seg;
	pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
	pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_plane(new pcl::PointCloud<pcl::PointXYZ>());
	seg.setOptimizeCoefficients(true);
	seg.setModelType(pcl::SACMODEL_PLANE);
	seg.setMethodType(pcl::SAC_RANSAC);
	seg.setMaxIterations(300);
	seg.setDistanceThreshold(0.02);

	int i = 0, nr_points = (int)cloud_filtered->size();
	int loop = 0;
	PointCloud<PointXYZ> cloud_plane_merged;
	//while (cloud_filtered->size() > 0.3 * nr_points)
	while(loop < 1)
	{
		// Segment the largest planar component from the remaining cloud
		seg.setInputCloud(cloud_filtered);
		seg.segment(*inliers, *coefficients);
		if (inliers->indices.size() == 0)
		{
			std::cout << "Could not estimate a planar model for the given dataset." << std::endl;
			break;
		}

		// Extract the planar inliers from the input cloud
		pcl::ExtractIndices<pcl::PointXYZ> extract;
		extract.setInputCloud(cloud_filtered);
		extract.setIndices(inliers);
		extract.setNegative(false);

		// Get the points associated with the planar surface
		extract.filter(*cloud_plane);
		cloud_plane_merged += *cloud_plane;
		std::cout << "PointCloud representing the planar component: " << cloud_plane->size() << " data points." << std::endl;

		// Remove the planar inliers, extract the rest
		extract.setNegative(true);
		extract.filter(*cloud_f);
		*cloud_filtered = *cloud_f;
		++loop;
	}
	std::cout << "planar point count:" << cloud_plane_merged.points.size() << std::endl;
	for (size_t i = 0; i < cloud_plane_merged.points.size(); ++i)
	{
		PointXYZRGB p;
		p.x = cloud_plane_merged.points[i].x;
		p.y = cloud_plane_merged.points[i].y;
		p.z = cloud_plane_merged.points[i].z;
		p.r = colors[0].r;
		p.g = colors[0].g;
		p.b = colors[0].b;
		cloud_out.points.push_back(p);
	}
	std::pair<int, int> plane;
	plane.first = 0;
	plane.second = cloud_plane_merged.points.size();
	clusters.push_back(plane);
	std::cout << "non planar point count:" << cloud_filtered->points.size() << std::endl;
	// Creating the KdTree object for the search method of the extraction
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(cloud_filtered);

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
	ec.setClusterTolerance(0.02); // 2cm
	ec.setMinClusterSize(100);
	ec.setMaxClusterSize(25000);
	ec.setSearchMethod(tree);
	ec.setInputCloud(cloud_filtered);
	ec.extract(cluster_indices);

	int j = 1;
	for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster(new pcl::PointCloud<pcl::PointXYZ>);
		for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); ++pit)
			cloud_cluster->push_back((*cloud_filtered)[*pit]); //*
		cloud_cluster->width = cloud_cluster->size();
		cloud_cluster->height = 1;
		cloud_cluster->is_dense = true;
		std::pair<int, int> rg;
		rg.first = clusters[j - 1].second;
		rg.second = rg.first + cloud_cluster->size();
		clusters.push_back(rg);
		for (int i = 0; i < cloud_cluster->points.size(); ++i)
		{
			PointXYZRGB p;
			p.x = cloud_cluster->points[i].x;
			p.y = cloud_cluster->points[i].y;
			p.z = cloud_cluster->points[i].z;
			p.r = colors[j].r;
			p.g = colors[j].g;
			p.b = colors[j].b;
			cloud_out.points.push_back(p);
		}
		if (j > 7)
			break;
		++j;
	}
}

void creatingClusters()
{
	cloud_clusters.push_back(cloud_total);
	for (int i = 0; i < clusters.size(); ++i)
	{
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
		for (size_t j = clusters[i].first; j < clusters[i].second; ++j)
		{
			cloud_ptr->points.push_back(cloud_total->points[j]);
		}
		cloud_clusters.push_back(cloud_ptr);
	}
}

void drawScene(Eigen::Affine3f &camPose, int x, int y, int width, int height, std::string &caption, boost::shared_ptr<pcl::PolygonMesh> mesh_ptr)
{
	_viewer = createViewer(camPose, x, y, width, height, caption);
	_mesh_ptr = mesh_ptr;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::fromPCLPointCloud2(_mesh_ptr->cloud, *cloud_in);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloud_out(new pcl::PointCloud<pcl::PointXYZRGB>);
	//std::vector<std::pair<int, int>> clusters;
	segmentPointCloud(cloud_in, *cloud_total, clusters);
	creatingClusters();
	while (loop)
	{
		_viewer->removePolygonMesh("mesh");
		_viewer->removeAllPointClouds();
		_viewer->removeText3D("length_picked");
		_viewer->addPointCloud<pcl::PointXYZRGB>(cloud_clusters[cluster_no], "clustering");
		//_viewer->addPointCloud<pcl::PointXYZRGB>(cloud_total, "clustering");
		//_viewer->addPolygonMesh(*mesh_ptr, "mesh");
		_viewer->spinOnce(10);
	}
	_viewer->close();
}