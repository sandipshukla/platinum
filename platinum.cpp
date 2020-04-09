#define _CRT_SECURE_NO_WORNINGS
#include <algorithm>
#include <condition_variable>
#include <atomic>
#include <cstring>
#include <experimental/filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <numeric>
#include <omp.h>
#include <bits/stdc++.h>
#include <optional>
#include <queue>
#include <random>
#include <thread>
#include <unordered_map>

using namespace std;

#ifdef _MSC_VER
namespace file_system = std::experimental::filesystem;
using integer = int;

//reverses byte order
integer byte_swap(integer val){return _byteswap_ulong(integer);}

// sanitize directory separator
string path(string path_string){return path_string;}

#elif defined(__GNUC__)
namespace file_system = std::experimental::filesystem;
using integer = int;

integer byte_swap(integer param_val){return __builtin_bswap32(param_val);}

string path(string param_path_string)
{
    replace(param_path_string.begin(), param_path_string.end(), '\\', '/');
    return param_path_string;
}
#endif

// aliases and constants
template <class T> using op = optional<T>;
using decimal = double;
using CHAR = char;
constexpr decimal INFINITE = 1e+10;
constexpr decimal EPSILON = 1e-4;
constexpr decimal PI = 3.14159265358979323846;

// vector maths
// 3d vector
template <typename T = decimal>
class vector3d
{
    protected:
        T _x;
        T _y;
        T _z;
    public:
        vector3d(T param_default = 0) : _x(param_default), _y(param_default), _z(param_default){}
        vector3d(T param_x, T param_y, T param_z) : _x(param_x), _y(param_y), _z(param_z){}
        ~vector3d(){}
        
        // get members
        const T x() const { return _x; }
        const T y() const { return _y; }
        const T z() const { return _z; }
    
        void set_x(T param_x) { _x = param_x; }
        void set_y(T param_y) { _y = param_y; }
        void set_z(T param_z) { _z = param_z; }
        
        vector3d<T> operator + (const vector3d<T>& param_vec) const { return vector3d<T>(_x + param_vec.x(), _y + param_vec.y(), _z + param_vec.z()); }
        vector3d<T> operator - (const vector3d<T>& param_vec) const { return vector3d<T>(_x - param_vec.x(), _y - param_vec.y(), _z - param_vec.z()); }
        vector3d<T> operator - () const { return vector3d<T>(-_x, -_y, -_z); }
        vector3d<T> operator * (const vector3d<T>& param_vec) const { return vector3d<T>(_x * param_vec._x, _y * param_vec._y, _z * param_vec._z); }
        vector3d<T> operator * (const decimal& param_scalar) const { return vector3d<T>(param_scalar * _x, param_scalar * _y , param_scalar * _z); }
        vector3d<T> operator / (const vector3d<T>& param_vec) const { return {_x / param_vec._x, _y / param_vec._y, _z / param_vec._z}; }
        vector3d<T> operator / (const T& param_scalar) const { return vector3d<T>(_x / param_scalar, _y / param_scalar, _z / param_scalar); }
        decimal operator [] (integer param_index) const {return (&_x)[param_index];}
        vector3d<T> operator = (const vector3d<T>& param_vec) 
        {   
            _x = param_vec._x;
            _y = param_vec._y;
            _z = param_vec._z;
            return *this;
        }
        
        
            
        T coord_max() const { return std::max({_x, _y, _z}); }
};

// Random number generator
class rng
{
    public:
        rng(){};
        rng(integer seed)
        {
            eng.seed(seed);
            uniform_real_dist.reset();
        }
        // Sample uniform random numbers in [0,1)
        decimal uniform_real(){ return uniform_real_dist(eng); }
        
        // Cosine-weighted direction sampling
        vector3d<decimal> cosine_weighted_direction()
        {
            decimal rand_sqrt = sqrt(uniform_real());
            decimal rotate_theta = 2 * PI * uniform_real();
            decimal rotate_x = rand_sqrt * cos(rotate_theta);
            decimal rotate_y = rand_sqrt * sin(rotate_theta);
            return vector3d<decimal>(rotate_x, rotate_y, sqrt(max(.0, 1 - rotate_x * rotate_x - rotate_y * rotate_y)));
        }
        
    private:
        mt19937 eng;
        uniform_real_distribution<decimal> uniform_real_dist;
};
        
// Square
template<typename T = decimal>
T square(T val){ return val * val;}

// Element wise min
template<typename T = vector3d<decimal>>
T vector_min( T param_vec1, T param_vec2) { return {min(param_vec1.x(), param_vec2.x()), min(param_vec1.y(), param_vec2.y()), min(param_vec1.z(), param_vec2.z()) }; }

// Element wise max
template<typename T = vector3d<decimal>>
T vector_max( T param_vec1, T param_vec2) { return {max(param_vec1.x(), param_vec2.x()), max(param_vec1.y(), param_vec2.y()), max(param_vec1.z(), param_vec2.z()) }; }

// Vector operations

// Dot product
template<typename T = vector3d<decimal>>
decimal dot(const T& param_vec1, const T& param_vec2){ return param_vec1.x() * param_vec2.x() + param_vec1.y() * param_vec2.y() + param_vec1.z() * param_vec2.z();}

// Cross product
template<typename T = vector3d<decimal>>
T cross(T param_vec1,T param_vec2)
{
    return
    {
        param_vec1.y() * param_vec2.z() - param_vec1.z() * param_vec2.y(),
        param_vec1.z() * param_vec2.x() - param_vec1.x() * param_vec2.z(),
        param_vec1.x() * param_vec2.y() - param_vec1.y() * param_vec2.x()
    };
}

// Reflected vector along normal vector
template<typename T = vector3d<decimal>>
T reflect(T param_wi, T param_normal){ return param_normal * dot(param_wi, param_normal) * 2 - param_wi; }


// Refracted vector along normal vector
template<typename T>
op<vector3d<decimal>> refract(T param_wi, T param_normal, decimal param_eta)
{
    const decimal var_t1 = dot(param_wi, param_normal);
    const decimal var_t2 = 1 - param_eta * param_eta * (1 - var_t1 * var_t1);
    return var_t2 > 0 ? ( param_normal * var_t1  - param_wi ) * param_eta -  param_normal * sqrt(var_t2) : op<vector3d<decimal>>{};
    
}

// Normalize
template<typename T = vector3d<>>
T normalize(T param_vec){ return param_vec / sqrt(dot(param_vec, param_vec));}

// Orthonormal basis compute [Duff et al. 2017]
template<typename T>
tuple<T, T> ortho_normal_basis(T param_normal)
{
    const decimal var_s = copysign(1, param_normal.z());
    const decimal var_a = -1 / (var_s + param_normal.z());
    const decimal var_b = param_normal.x() * param_normal.y() * var_a;
    const T vec_u(1 + var_s * param_normal.x() * param_normal.x() * var_a, var_s * var_b, -var_s * param_normal.x());
    const T vec_v(var_b, var_s + param_normal.y() * param_normal.y() * var_a, -param_normal.y());
    return {vec_u,vec_v};
}

// rays
template <typename T = vector3d<>>
struct ray
{
    T origin;
    T direction;
};

// Axis aligned bounding box
template<typename T = vector3d<>>
class aabb
{
        
    private:
        T _minimum = T(INFINITE);
        T _maximum = T(-INFINITE);
    public:
        aabb(T param_minimum, T param_maximum) : _minimum(param_minimum), _maximum(param_maximum){}
        aabb() {}
        const T operator[](integer param_int) const { return (&_minimum)[param_int];}
                
        const T centroid() const { return (_minimum + _maximum) * 0.5; }
        const T minimum() const { return _minimum; }
        const T maximum() const { return _maximum; }
                
        decimal surface_area() const
        {
            T length_vector = _maximum - _minimum;
            return 2 * (length_vector.x() * length_vector.y() + length_vector.y() * length_vector.z() + length_vector.z() * length_vector.x());
        }
        
        // test aabb intersection with ray from <http://psgraphics.blogspot.com/2016/02/new-simple-ray-box-test-from-andrew.html>
        bool intersect_ray(const ray<T>& param_ray, decimal param_tmin, decimal param_tmax) const
        {         
            for (integer index = 0; index < 3; index++)
            {
                const decimal inverse_direction = 1 / param_ray.direction[index];
                decimal var_tmin = (_minimum[index] - param_ray.origin[index]) * inverse_direction;
                decimal var_tmax = (_maximum[index] - param_ray.origin[index]) * inverse_direction;
                if( inverse_direction < 0){ swap(var_tmin, var_tmax); }
                param_tmin  = max(var_tmin, param_tmin);
                param_tmax  = min(var_tmax, param_tmax);
                if( param_tmax < param_tmin )
                {
                    return false;                    
                }
            }
            return true;
        }
};

// merge aabb with a point
template <typename T = vector3d<>>
aabb<T> merge_aabb( aabb<T> param_aabb, T param_point ){ return { vector_min<T>(param_aabb.minimum(), param_point), vector_max<T>(param_aabb.maximum(), param_point)}; }
// merge two aabb
template <typename T = vector3d<>>
aabb<T> merge_aabb( aabb<T> param_aabb1, aabb<T> param_aabb2 ){ return {vector_min<T>(param_aabb1.minimum(), param_aabb2.minimum()), vector_max<T>(param_aabb1.maximum(), param_aabb2.maximum())}; }
    
// 1d discrete distribution
template <typename T = decimal>
class dist1d
{
public:
    
    const vector<T>& cdf() const { return _cdf; }
    // Add a value to the distribution
    void add(T param_val){ _cdf.push_back(_cdf.back() + param_val); }
    
    // normalize the distribution
    void normalize()
    {
        auto sum = _cdf.back();
        for ( T& val : _cdf)
        {
            val /= sum;
        }
        
    }
    
    // evaluate pmf
    T pmf(integer param_val) const
    {
        return ( param_val < 0 || param_val + 1 >= integer(_cdf.size())) ? 0 : _cdf[param_val + 1] - _cdf[param_val];
    }
    
    // sample the distribution
    integer sample(rng& param_rng) const
    {
        decimal var_rng_num = param_rng.uniform_real();
        const auto iter = upper_bound(_cdf.begin(), _cdf.end(), var_rng_num);
        return clamp(integer(distance(_cdf.begin(), iter)) - 1, 0, integer(_cdf.size()) - 2);
    }
    
private:
    vector<T> _cdf{0};
};
    
// 2d discrete distribution
template<typename T = decimal>
class dist2d
{
    public:
        // add values to distribution
        void init(const vector<T>& param_vals, integer param_width, integer param_height)
        {
            width = param_width;
            height = param_height;
            conditional_dist.assign(height, {});
            for (integer index1 = 0; index1 < height; index1++)
            {
                auto& dist = conditional_dist[index1];
                for (integer index2 = 0; index2 < width; index2++)
                {
                    dist.add(param_vals[index1 * width + index2]);
                }
                marginal_dist.add(dist.cdf().back());
                dist.normalize();
            }
            marginal_dist.normalize();
        }
        
        // Evaluate pmf
        T pmf(T param_u, T param_v) const
        {
            const integer index_height = min(integer(param_u * param_v), height - 1);
            return marginal_dist.pmf(index_height) * conditional_dist[index_height].pmf(integer(param_u * width)) * width * height;
        }
        
        // sample the distribution
        tuple<T, T> sample(rng& param_rng) const
        {
            const integer index_y = marginal_dist.sample(param_rng);
            const integer index_x = conditional_dist[index_y].sample(param_rng);
            return { (index_x + param_rng.uniform_real()) / width, (index_y + param_rng.uniform_real()) / height};
        }
    private:
        // conditional distribution corresponding to a row
        vector<dist1d<T>> conditional_dist;
        //Marginal distribution
        dist1d<T> marginal_dist;
        // size of the distribution
        integer width, height;
    
};

template <typename T = vector3d<decimal>>
T interpolate_barycentric(T param_vertex_a, T param_vertex_b, T param_vertex_c, decimal param_interpolator_u, decimal param_interpolator_v)
{
    return param_vertex_a * (1 - param_interpolator_u - param_interpolator_v) + param_vertex_b * param_interpolator_u + param_vertex_c * param_interpolator_v;
}

struct hit_point_triangle
{
    decimal distance;
    decimal barycentric_u;
    decimal barycentric_v;
};

template <typename T = vector3d<decimal>>
class triangle
{
    private:
        T _vertex;
        T _edge1;
        T _edge2;
        T _normal;
        aabb<T> _bound;
        T _bound_centroid;
        integer _object_index;
        integer _face_index;
        integer id;
    public:
        triangle( integer param_id, const T param_vertex1, const T param_vertex2, const T param_vertex3, integer param_object_index, integer param_face_index) : id(param_id), _vertex(param_vertex1), _object_index(param_object_index), _face_index(param_face_index)
        {
            _edge1 = param_vertex2 - param_vertex1;
            _edge2 = param_vertex3 - param_vertex1;
            _normal = normalize(cross(_edge1, _edge2));
            _bound = merge_aabb<>(_bound, param_vertex1);
            _bound = merge_aabb<>(_bound, param_vertex2);
            _bound = merge_aabb<>(_bound, param_vertex3);
            _bound_centroid = _bound.centroid();
        }
        
        const T& normal() const { return _normal; }
        const aabb<T>& bound() const { return _bound; }
        const T& bound_centroid() const { return _bound_centroid; }
        const integer object_index() const { return _object_index; }
        const integer face_index() const { return _face_index; }
    
        //Moller & Trumbore 1997
        op<hit_point_triangle> intersect_ray(const ray<T>& param_ray, decimal param_distance_low, decimal param_distance_high) const
        {
            T var_re2 = cross(param_ray.direction, _edge2);
            T var_rv = param_ray.origin - _vertex;
            T var_cross_rv_e1 = cross(var_rv, _edge1);
            decimal var_dot_e1_re2 = dot(_edge1, var_re2);
            decimal var_abs_dot_e1_re2 = abs(var_dot_e1_re2);
            decimal var_dot_sign = copysign(1, var_dot_e1_re2);
            decimal var_u = dot(var_rv, var_re2) * var_dot_sign;
            decimal var_v = dot(param_ray.direction, var_cross_rv_e1) * var_dot_sign;
            if (var_abs_dot_e1_re2 < 1e-8 || var_u < 0 || var_v < 0 || var_u + var_v > var_abs_dot_e1_re2) { return{}; }
            decimal var_distance = dot(_edge2, var_cross_rv_e1) / var_dot_e1_re2;
            return
                var_distance < param_distance_low || param_distance_high < var_distance ? op<hit_point_triangle>{} : hit_point_triangle{var_distance, var_u / var_abs_dot_e1_re2, var_v / var_abs_dot_e1_re2};
        }

        op<hit_point_triangle> isect(const ray<T>& r, decimal tl, decimal th) const 
        {
            vector3d<> p = cross(r.direction, _edge2), tv = r.origin - _vertex;
            vector3d<> q = cross(tv, _edge1);
            decimal d = dot(_edge1, p), ad = abs(d), s = copysign(1, d);
            decimal u = dot(tv, p) * s, v = dot(r.direction, q) * s;
            if (ad < 1e-8 || u < 0 || v < 0 || u + v > ad) {
                return {};
            }
            decimal t = dot(_edge2, q) / d;
            return t < tl || th < t ? op<hit_point_triangle>{} : hit_point_triangle{t, u / ad, v / ad};
        }
};

struct bvh_node
{
    bool leaf = false;
    aabb<> bound;
    integer triangle_indices_start;
    integer triangle_indices_end;
    integer index_child_left;
    integer index_child_right;
};

// surface geometry data
template <typename T = vector3d<decimal>>
struct geometry
{
    vector<T> positions;
    vector<T> normals;
    vector<T> tex_coords;
};

// indices for face data
struct face_index
{
    integer position = -1;
    integer normal = -1;
    integer tex_coord = -1;
};
    
//surface point
template <typename T = vector3d<decimal>>
class surface_point
{
    public:
        surface_point(){}
        surface_point(T param_position, T param_normal, T param_tex_coord) : _position(param_position), _normal(param_normal), _tex_coord(param_tex_coord)
        {
            tie(_tangent_u, _tangent_v) = ortho_normal_basis(_normal);
        }
        
        // returns true if wi and wo are in the same direction according to the normal
        bool op(T param_wi, T param_wo) const
        {
            return dot(param_wi, _normal) * dot(param_wo, _normal) <= 0;
        }
        
        T position() const { return _position; }
        T position() { return _position; }
        T normal() const { return _normal; }
        T normal() { return _normal; }
        T tex_coord()  const { return _tex_coord; }
        T tex_coord()  { return _tex_coord; }
        T tangent_u() const { return _tangent_u; }
        T tangent_u() { return _tangent_u; }
        T tangent_v() const { return _tangent_v; }
        T tangent_v() { return _tangent_v; }
        
        // returns orthonormal basis according to the incident direction wi at the surface point
        tuple<T, T, T> ortho_normal_basis_wi(T param_wi) const
        {
            const integer var_incident_factor = dot(param_wi, _normal) > 0;
            return { var_incident_factor ? _normal : -_normal , _tangent_u, var_incident_factor ? _tangent_v : -_tangent_v };
        }
        
    private:
        T _position;
        T _normal;
        T _tex_coord;
        T _tangent_u, _tangent_v;
};

// hit point
struct hit_point
{
    surface_point<vector3d<decimal>> surface_hit_point;
    struct object* object_hit_ptr;
};

// geometry factor between surfaces
decimal geometric_factor(const surface_point<vector3d<decimal>>& param_surface_point1, const surface_point<vector3d<decimal>>& param_surface_point2)
{
    auto var_direction = param_surface_point2.position() - param_surface_point1.position();
    const decimal var_l2 = dot(var_direction, var_direction);
    var_direction = var_direction / sqrt(var_l2);
    return abs(dot(param_surface_point1.normal(), var_direction)) * abs(dot(param_surface_point2.normal(), -var_direction)) / var_l2;
}

// texture
class texture
{
    public:
        texture(): _width(0), _height(0){}
        void name(const string param_name) { _name = param_name; }
        const vector<decimal>& colors() const { return _colors; }
        const vector<decimal>& alphas() const { return _alphas; }
        const decimal get_width() const { return _width; }
        const decimal get_height() const { return _height; }

        //vertically flipped pixel coordinate
        integer flipped_coordinate(integer param_index)
        {
            const integer var_j = param_index / 3;
            const integer var_x = var_j % _width;
            const integer var_y = var_j / _width;
            return 3 * ( (_height - var_y - 1) * _width + var_x ) + param_index % 3;
        }
        
        // post process pixel for pmp textures
        decimal post_process(integer param_index, decimal param_exp, vector<uint8_t>& param_texture_colors)
        {
            decimal var_color = decimal(param_texture_colors[flipped_coordinate(param_index)]);
            return pow( var_color / param_exp, 2.2);
        }
        
        // process pixel for pmf textures
        decimal post_process(integer param_index, decimal param_exp, vector<decimal>& param_texture_colors)
        {
            if (param_exp < 0)
            {
                return param_texture_colors[flipped_coordinate(param_index)];
            }
            auto var_m = byte_swap(*(int32_t *)&param_texture_colors[flipped_coordinate(param_index)]);
            return *(float*)&var_m;
        }

        // load textures.
        template<class T> void load(vector<decimal>& param_texture_colors, string param_file_name)
        {
            static vector<T> var_texture_colors;
            const char* var_file_tex = param_file_name.c_str();
            FILE* tex_file = fopen(param_file_name.c_str(), "rb");
            if (!tex_file){ return; }
            
            decimal var_exp;
            char* var_ppm_format;
            fscanf(tex_file, "%*s %d %d %lf%*c", &_width, &_height, &var_exp);
            const integer tex_size = _width * _height * 3;
            var_texture_colors.assign(tex_size, 0);
            _colors.assign(tex_size, 0);
            fread(var_texture_colors.data(), sizeof(T), tex_size, tex_file);
            for (integer tex_index = 0; tex_index < tex_size; tex_index++)
            {
                _colors[tex_index] = post_process(tex_index, var_exp, var_texture_colors);
            }
            fclose(tex_file);
        }
        
        // load pfm texture
        void loadpfm(string param_file_name){ load<decimal>(_colors, param_file_name); }
        // load ppm texture
        void loadppm(string param_file_name)
        {
            auto var_file_path = file_system::path(param_file_name);
            auto var_file_path_color = var_file_path.replace_extension(".ppm").string();
            auto var_file_path_alpha = var_file_path.replace_extension("_alpha.ppm").string();
            load<uint8_t>(_colors, var_file_path_color);
            load<uint8_t>(_alphas, var_file_path_alpha);
        }
        
        // evaluate texture at the given pixel coordinate
        vector3d<decimal> eval(vector3d<decimal> param_tex_coord, bool alpha = 0) const
        {
            const decimal var_tex_coord_u = param_tex_coord.x() - floor(param_tex_coord.x());
            const decimal var_tex_coord_v = param_tex_coord.y() - floor(param_tex_coord.y());
            const integer var_width = clamp(integer(var_tex_coord_u * _width), 0, _width - 1);
            const integer var_height = clamp(integer(var_tex_coord_v * _height), 0, _height - 1);
            const integer var_index = _width * var_height + var_width;
            return alpha ? vector3d<decimal>(_alphas[3 * var_index]) : vector3d<decimal>(_colors[3 * var_index], _colors[3 * var_index + 1], _colors[3 * var_index + 2]);
        }
        
    private:
        string _name;
        integer _width;
        integer _height;
        vector<decimal> _colors;
        vector<decimal> _alphas;
};

// Lights and material types
namespace lnm
{
    enum
    {
        none = 0,    // Uninitialized
        area_light = 1 << 0,
        environment_light = 1 << 1,
        sensor = 1 << 2,
        diffuse_material = 1 << 3,
        glossy_material = 1 << 4,
        transparent_mask = 1 << 5,
        fresnel_reflection = 1 << 6,
        fresnel_transmission = 1 << 7,
        mirror_reflection = 1 << 8,
        
        // Compound Types
        transmissive_material = fresnel_transmission | transparent_mask,
        fresnel = fresnel_reflection | fresnel_transmission,
        lights = area_light | environment_light,
        specular_material = mirror_reflection | fresnel | transparent_mask,
        non_specular_material = diffuse_material | glossy_material,
        reflection = diffuse_material | glossy_material | mirror_reflection | fresnel_reflection,
        
    };
    
}

template<typename T = decimal>
class material
{
    private:
        integer _type = lnm::none;
        vector3d<T> _kd;       //diffuse reflectance
        vector3d<T> _ks;       // specular reflectance
        vector3d<T> _ke;       // luminance
        texture* _map_kd = 0;
        T _ior;                //index of refraction
        T _specular_exponent;       // for phong shading
        T _anisotropy;
        T _roughness_x;
        T _roughness_y;
    
    public:
        integer type_ = lnm::none;
        material(integer param_lnm_type = lnm::none, vector3d<> param_kd = vector3d<>(1)) : _type(param_lnm_type), _kd(param_kd){}
    
        void type(integer param_lnm_type ) { _type = param_lnm_type; }
        void kd(const vector3d<> param_kd) { _kd = param_kd; }
        void ks(const vector3d<> param_ks) { _ks = param_ks; }
        void ke(const vector3d<> param_ke) { _ke = param_ke; }
        void map_kd(texture* param_ptr_map) { _map_kd = param_ptr_map; }
        void ior(const decimal param_ior) { _ior = param_ior; }
        void specular_exponent(const decimal param_specular_exponent) { _specular_exponent = param_specular_exponent; }
        void anisotropy(const decimal param_anisotropy) { _anisotropy = param_anisotropy; }
        void roughness_x(const decimal param_roughness_x) { _roughness_x = param_roughness_x; }
        void roughness_y(const decimal param_roughness_y) { _roughness_y = param_roughness_y; }
    
        const integer type() const { return _type; }
        const vector3d<T> kd() const { return _kd; }
        const vector3d<T> ks() const { return _ks; }
        const vector3d<T> ke() const { return _ke; }
        const texture* map_kd() const { return _map_kd; }
        const T ior() const { return _ior; }
        const T specular_exponent() const { return _specular_exponent; }
        const T anisotropy() const { return _anisotropy; }
        const T roughness_x() const { return _roughness_x; }
        const T roughness_y() const { return _roughness_y; }
        
        const T pdf(integer param_object_type, const surface_point<>& param_surface_point, const vector3d<decimal>& param_wi, const vector3d<decimal>& param_wo) const
        {
            if (param_surface_point.op(param_wi, param_wo)) { return 0; }
            if (param_object_type & lnm::diffuse_material) { return 1 / PI; }
            else if (param_object_type & lnm::glossy_material)
            {
                vector3d<decimal> var_wh = normalize(param_wi + param_wo);
                const auto[param_normal, param_tangent, param_bi_tangent] = param_surface_point.ortho_normal_basis_wi(param_wi);
                return aniso_ggx_distribution(var_wh, param_tangent, param_bi_tangent, param_normal) * dot(var_wh, param_normal) / (4 * dot(param_wo, var_wh) * dot(param_wo, param_normal));
            }
            return 0;
            
        }
        
        // Anisotropic  GGX Normal Distributions(Burley 2012)
        T aniso_ggx_distribution(vector3d<T> param_half_vector, vector3d<T> param_vector_u, vector3d<T> param_vector_v, vector3d<T> param_normal) const
        {
            return 1 / (PI*_roughness_x*_roughness_y*square(square(dot(param_half_vector, param_vector_u)/_roughness_x) + square(dot(param_half_vector, param_vector_v)/_roughness_y) + square(dot(param_half_vector, param_normal))));
        }
        // Smith's Geometric term for anisotropic GGX
        T smith_g_aniso_ggx(vector3d<decimal> param_wi, vector3d<decimal> param_wo, vector3d<decimal> param_u, vector3d<decimal> param_v, vector3d<decimal> param_normal) const
        {
            auto smith_g1 = [&](vector3d<decimal> param_w)
            {
                const decimal var_w_dot_n = dot(param_w, param_normal);
                const decimal var_sqrt = sqrt(1 - var_w_dot_n);
                const decimal var_w_dot_u = dot(param_w, param_u) / var_sqrt;
                const decimal var_w_dot_v = dot(param_w, param_v) / var_sqrt;
                const decimal var_a_squared = square(var_w_dot_u * _roughness_x) + square(var_w_dot_v * _roughness_y);
                return var_w_dot_n == 0 ? 0 : 2 / (1 + sqrt(1 + var_a_squared * square( var_sqrt / var_w_dot_n)));
            };
            return smith_g1(param_wi) * smith_g1(param_wo);
        }
               
};

// direction sampling
struct direction_sample
{
    integer object_type;
    ray<vector3d<decimal>> ray_sampled;
    vector3d<decimal> evaluated_w;
    decimal prob_component_select;
};

// light sampling result
struct light_sample
{
    vector3d<decimal> wo;   // sampled direction
    decimal sample_distance;
    vector3d<decimal> evaluated_le;
    decimal evaluated_prob;
};

// Kolb et al. 1995
struct lens
{
    decimal curvature_radius;
    decimal thickness;
    decimal eta;
    decimal aperture_radius;
};

struct lens_hit_info
{
    decimal distance; // distance along the ray
    vector3d<decimal> normal; // hit normal
};

struct sensor
{
    vector3d<decimal> position;
    vector3d<> basis_u;
    vector3d<> basis_v;
    vector3d<> basis_w;
    tuple< vector3d<decimal>, vector3d<decimal>, vector3d<decimal> > basis_uvw;
    decimal focus_distance;
    decimal aspect_ratio;
    decimal diagonal_length;
    decimal distance_nearest_element;
    decimal sensitivity;
    vector<lens> lens_elements;
    vector<op<aabb<vector3d<decimal>>>> bounds_exit_pupil;       // for importance sampling
    
};

struct light_sampling_params
{
    dist1d<decimal> discrete_dist1d;
    decimal inverse_area;
};

struct env_light_sampling_params
{
    texture map;
    dist2d<decimal> discrete_dist2d;
    decimal rotation;                       // rotation of map around (0, 1, 0)
};

class object
{
    private:
        integer _type;
        vector<face_index> _face_indices;
        material<decimal>* _material = 0;
        sensor _sensor;
        light_sampling_params _light_sampling_params;
        env_light_sampling_params _env_light_sampling_params;

    public:        
        object(integer param_lnm_type = lnm::none, material<>* param_ptr_material = nullptr) : _type(param_lnm_type), _material(param_ptr_material){}
        const integer type() const { return _type; }
        void type(const integer& param_type) { _type = param_type; }
        vector<face_index>& face_indices() { return _face_indices;}
        material<decimal>* get_material() { return _material;}

        // realistic lens system from PBRT->6.4(3rd edition)
        op<ray<vector3d<decimal>>> trace_lens_system(ray<vector3d<decimal>> param_ray) const
        {
            decimal var_z  = 0;
            for( auto index = integer(_sensor.lens_elements.size()) - 1; index >= 0; index--)
            {
                const auto& var_lens_element = _sensor.lens_elements[index];
                var_z -= var_lens_element.thickness;
            
                auto var_lens_hit_info = [&]() -> op<lens_hit_info>
                {
                    // aperture stop
                    if (var_lens_element.curvature_radius == 0)
                    {
                        // check intersection with aperture stop
                        decimal var_t = (var_z - param_ray.origin.z()) / param_ray.direction.z();
                        return var_t < 0 ? op<lens_hit_info>{} : lens_hit_info{var_t, {}};
                    }
                    
                    // case spherical lens element intersection
                    vector3d<decimal> var_z_center(0, 0, var_z + var_lens_element.curvature_radius);
                    const vector3d<decimal> var_object_center = var_z_center - param_ray.origin;
                    const decimal var_intersect_factor = dot<vector3d<decimal>>(var_object_center, param_ray.direction);
                    const decimal var_dt = var_intersect_factor * var_intersect_factor - dot<vector3d<decimal>>(var_object_center, var_object_center) + var_lens_element.curvature_radius * var_lens_element.curvature_radius;
                    if( var_dt < 0)
                    {
                        // no intersection
                        return {};
                        
                    }
                    
                    const decimal var_t0 = var_intersect_factor - sqrt(var_dt);
                    const decimal var_t1 = var_intersect_factor + sqrt(var_dt);
                    const decimal var_t = (param_ray.direction.z() > 0) ^ ( var_lens_element.curvature_radius < 0) ? min(var_t0, var_t1) : max(var_t0, var_t1);
                    
                    if (var_t < 0)
                    {
                        // no intersection towards positive direction
                        return {};
                    }
                    
                    vector3d<decimal> var_normal = ( param_ray.origin +  param_ray.direction * var_t - var_z_center) / var_lens_element.curvature_radius;
                    var_normal = dot<vector3d<decimal>>(var_normal, -param_ray.direction);
                    return lens_hit_info{var_t, var_normal};
                    
                }();
                if (!var_lens_hit_info){ return {}; }
                
                // intersection with aperture
                const vector3d<decimal> hit_point = param_ray.origin + param_ray.direction * var_lens_hit_info->distance;
                if( hit_point.x() * hit_point.x() + hit_point.y() * hit_point.y() > var_lens_element.aperture_radius * var_lens_element.aperture_radius)
                {
                    return {};
                }
                
                // next ray setup
                // aperture stop
                param_ray.origin = hit_point;
                if( var_lens_element.curvature_radius == 0)
                {
                    // continue along the direcion
                    continue;
                }
                
                // refraction by spherical lens element
                const decimal var_eta_i = var_lens_element.eta / ( index > 0 && _sensor.lens_elements[index-1].eta != 0 ? _sensor.lens_elements[index-1].eta : 1);
                const auto var_eta_t = refract<vector3d<decimal>>(-param_ray.direction, var_lens_hit_info->normal, var_eta_i);
                if(!var_eta_t)
                {
                    // total internal reflection
                    return{};
                }
                param_ray.direction = *var_eta_t;
            }
            return param_ray;
        }
    
        // calculate focus distance based on distance between the last lens element and the image plane
        decimal focus_distance( decimal param_intersect_distance) const
        {
            op<ray<vector3d<decimal>>> var_ray;
            const auto& var_lens_back = _sensor.lens_elements.back();
            // shoot rays parellel to the optical axis
            for (integer index = 9; index > 0; index--)
            {
                if ( var_ray = trace_lens_system({vector3d<decimal>(0, 0, -var_lens_back.thickness + param_intersect_distance), normalize<vector3d<decimal>>(vector3d<decimal>(var_lens_back.aperture_radius*index/10, 0, -param_intersect_distance))})){ break; }
            }
            
            if(!var_ray){ return INFINITE;}
            const decimal var_t = -var_ray->origin.x() / var_ray->direction.x();
            const decimal var_z = (var_ray->origin + var_ray->direction * var_t).z();
            decimal var_sz = 0;
            for( auto& var_lens_element : _sensor.lens_elements) { var_sz += var_lens_element.thickness; }
            
            return var_z < var_sz ? -var_z : INFINITE;
        }
    
        // sensor object initialization
        void init_sensor(string lens_file, vector3d<decimal> param_position, vector3d<decimal>  param_center, vector3d<decimal> param_u, decimal param_fov, decimal param_focus_distance, decimal param_diagonal_length, decimal param_sensitivity, decimal param_aspect_ratio)
        {
            cout << "\n" << "sensor initialization started" << "\n";
            _sensor.position = param_position;
            _sensor.focus_distance = tan(param_fov * PI / 180 * 0.5);
            _sensor.aspect_ratio = param_aspect_ratio;
            _sensor.diagonal_length = param_diagonal_length * .001;
            _sensor.sensitivity = param_sensitivity;
            vector3d<> var_w = normalize<vector3d<decimal>>(param_position - param_center);
            vector3d<> var_u = normalize<vector3d<decimal>>(cross<vector3d<decimal>>(param_u, var_w));
            vector3d<> var_v = cross<vector3d<decimal>>(var_w, var_u);
            _sensor.basis_u = var_u;
            _sensor.basis_v = var_v;
            _sensor.basis_w = var_w;
            
            // load lens data
            CHAR lens_data[4096];
            ifstream stream_lens_file(lens_file);
            while (stream_lens_file.getline(lens_data, 4096))
            {
                if (lens_data[0] == '#' || lens_data[0] == '\0'){ continue; }
                decimal var_curvature_radius, var_thickness, var_eta, var_aperture_radius;
                sscanf(lens_data, "%1f %1f %1f %1f", &var_curvature_radius, var_thickness, &var_eta, &var_aperture_radius);
                _sensor.lens_elements.push_back({var_curvature_radius * .001, var_thickness * .001, var_eta * .001, var_aperture_radius * .001 * .5 });                    
            }
            if (_sensor.lens_elements.empty()){ return; }
            
            // auto focus region
            decimal var_focus_distance_low = EPSILON, var_focus_distance_high = INFINITE;
            for (integer index = 0; index < 99; index++)
            {
                decimal var_auto_focus_distance = (var_focus_distance_low + var_focus_distance_high) * 0.5;
                (focus_distance(var_auto_focus_distance) < param_focus_distance ? var_focus_distance_high : var_focus_distance_low) = var_auto_focus_distance;
            }
            _sensor.distance_nearest_element = var_focus_distance_high;
            
            // compute exit pupil bounds at sampled positions on the image plane
            rng var_rng(137);
            integer var_num_bounds_exit_pupil = 64;
            _sensor.bounds_exit_pupil.assign(var_num_bounds_exit_pupil, {});
            decimal var_center_fov = -1;
            decimal var_sensor_y = sqrt(_sensor.diagonal_length * _sensor.diagonal_length / (1 + _sensor.aspect_ratio)), var_sensor_x = _sensor.aspect_ratio * var_sensor_y;
            integer var_m = 1 << 12;
            for (integer index1 = 0; index1 < var_num_bounds_exit_pupil; index1++)
            {
                aabb<vector3d<decimal>> var_bound_exit_pupil;
                bool var_lens_traced = false;
                auto& var_lens_element_back = _sensor.lens_elements.back();
                vector3d<decimal> var_distance_element_lens_last(0, 0, -var_lens_element_back.thickness + _sensor.distance_nearest_element);
                for (integer index2 = 0; index2 < var_num_bounds_exit_pupil; index2++)
                {
                    var_distance_element_lens_last.set_x((index1 + var_rng.uniform_real()) / var_num_bounds_exit_pupil * _sensor.diagonal_length * .5);
                    const decimal var_rng_sqrt = sqrt(var_rng.uniform_real());
                    const decimal var_angle = 2 * PI * var_rng_sqrt;
                    const vector3d<decimal> var_position_lens(var_rng_sqrt * cos(var_angle) * var_lens_element_back.aperture_radius, var_rng_sqrt * sin(var_angle) * var_lens_element_back.aperture_radius, -var_lens_element_back.thickness);
                    const auto& var_traced_ray = trace_lens_system({var_distance_element_lens_last, normalize<vector3d<decimal>>(var_position_lens - var_position_lens)});
                    if (var_traced_ray)
                    {
                        var_lens_traced = true;
                        var_bound_exit_pupil = merge_aabb<vector3d<decimal>>(var_bound_exit_pupil, var_position_lens);
                        if (var_distance_element_lens_last.x() < var_sensor_x * .5){ var_center_fov = max(var_center_fov, var_traced_ray->direction.z()); }
                    }
                }
                if (var_lens_traced) { _sensor.bounds_exit_pupil[index1] = var_bound_exit_pupil; }
            }
            // effective vertical feld of view
            //cout << "\n" (atan(tan(PI - acos(var_center_fov)) / _sensor.aspect_ratio) * 180 / PI * 2);
            
        }

        // setup objct as environment light
        void init_env_light(string param_env_map_file_name, decimal param_rotation)
        {
            _env_light_sampling_params.map.name(param_env_map_file_name);
            _env_light_sampling_params.map.loadpfm(param_env_map_file_name);
            _env_light_sampling_params.rotation = param_rotation;
            auto& var_env_map_colors = _env_light_sampling_params.map.colors();
            integer var_map_width = _env_light_sampling_params.map.get_width();
            integer var_map_height = _env_light_sampling_params.map.get_height();
            integer var_width_height_product = var_map_width * var_map_height;
            vector<decimal> var_light_sampling_params(var_width_height_product);
            
            for (integer index = 0; index < var_width_height_product; index++)
            {
                vector3d<decimal> var_color(var_env_map_colors[3 * index], var_env_map_colors[3 * index + 1], var_env_map_colors[3 * index + 2]);
                var_light_sampling_params[index] = var_color.coord_max() * sin(PI * (index / var_map_width + 0.5) / var_map_height);
            }
            _env_light_sampling_params.discrete_dist2d.init(var_light_sampling_params, var_map_width, var_map_height);
        }
    
        // setup object as area light
        void init_area_light(const geometry<vector3d<decimal>>& param_geometry)
        {
            for (size_t face_index = 0; face_index < _face_indices.size(); face_index += 3)
            {
                const vector3d<decimal> var_vertex1 = param_geometry.positions[_face_indices[face_index].position];
                const vector3d<decimal> var_vertex2 = param_geometry.positions[_face_indices[face_index + 1].position];
                const vector3d<decimal> var_vertex3 = param_geometry.positions[_face_indices[face_index + 2].position];
                const vector3d<decimal> var_vertex_normal = cross(var_vertex2 - var_vertex1, var_vertex3 - var_vertex1);
                _light_sampling_params.discrete_dist1d.add(sqrt(dot(var_vertex_normal, var_vertex_normal)) * 0.5);                
            }
            _light_sampling_params.inverse_area = 1 / _light_sampling_params.discrete_dist1d.cdf().back();
            _light_sampling_params.discrete_dist1d.normalize();
        }
    
        // directional sampling
        op<direction_sample> sample_direction(rng& param_rng, integer param_lnm_type, const surface_point<vector3d<decimal>>& param_surface_point, const vector3d<decimal>& param_wi)
        {
            if (param_lnm_type & lnm::sensor)
            {
                const vector3d<decimal> var_ray_point = param_wi * 2 - 1;
                const integer var_total_bounds_exit_pupil = integer(_sensor.bounds_exit_pupil.size());
                const auto& var_lens_back = _sensor.lens_elements.back();
                
                // pin hole camera
                if (_sensor.lens_elements.empty())
                {
                    // sensor co-ordinates directions
                    const vector3d<decimal> var_direction = -normalize(vector3d<decimal>(var_ray_point.x() * _sensor.aspect_ratio * _sensor.focus_distance, var_ray_point.y() * _sensor.focus_distance, 1));
                    return direction_sample{ lnm::sensor, _sensor.position, _sensor.basis_u * var_direction.x() + _sensor.basis_v * var_direction.y() + _sensor.basis_w * var_direction.z(), vector3d<decimal>(1), 1};
                }
                
                // realistic camera
                
                // find a position on sensor plane
                const decimal var_sensor_val = sqrt(_sensor.diagonal_length * _sensor.diagonal_length / (1 + _sensor.aspect_ratio * _sensor.aspect_ratio));
                const vector3d<decimal> var_sensor_origin = var_ray_point * vector3d<decimal>(_sensor.aspect_ratio * var_sensor_val * .5, var_sensor_val * .5, 0) + vector3d<decimal>(0, 0, _sensor.distance_nearest_element);
                
                // find bound corresponding to pixel location
                const decimal var_sqrt_sensor_origin = sqrt(var_sensor_origin.x() * var_sensor_origin.x() + var_sensor_origin.y() * var_sensor_origin.y());
                const integer var_index = clamp(integer(var_sqrt_sensor_origin / _sensor.diagonal_length * 2 * var_total_bounds_exit_pupil), 0, var_total_bounds_exit_pupil - 1);
                const auto& var_bound_exit_pupil = _sensor.bounds_exit_pupil[var_index];
                if(!var_bound_exit_pupil){ return{}; }
                
                // sample a location on exit pupil and generate initial ray direction
                const vector3d<decimal> var_bound_diagonal_vector = var_bound_exit_pupil->maximum() - var_bound_exit_pupil->minimum();
                const vector3d<decimal> var_position_exit_pupil = var_bound_exit_pupil->minimum() + var_bound_diagonal_vector * vector3d<decimal>(param_rng.uniform_real() , param_rng.uniform_real(), 0);
                const decimal var_sensor_plane_x = var_sensor_val != 0 ? var_sensor_origin.x() / var_sqrt_sensor_origin : 1;
                const decimal var_sensor_plane_y = var_sensor_val != 0 ? var_sensor_origin.y() / var_sqrt_sensor_origin : 0;
                const vector3d<decimal> var_direction = normalize(vector3d<decimal>(var_sensor_plane_x * var_position_exit_pupil.x() - var_sensor_plane_y * var_position_exit_pupil.y(), var_sensor_plane_y * var_position_exit_pupil.x() + var_sensor_plane_x * var_position_exit_pupil.y(), var_position_exit_pupil.z()) - var_sensor_origin);
                
                // now trace the rays through the lens
                const auto var_ray = trace_lens_system({var_sensor_origin, var_direction});
                if(!var_ray){ return {}; }
                
                // evaluate contribution
                const decimal var_bound_xy = var_bound_diagonal_vector.x() * var_bound_diagonal_vector.y();
                const decimal var_thickness = var_lens_back.thickness + _sensor.distance_nearest_element;
                const decimal var_factor_direction_z = var_direction.z() * var_direction.z() * var_direction.z() * var_direction.z() * var_bound_xy / var_thickness * var_thickness;
                
                return
                direction_sample
                {
                    lnm::sensor,
                    get<0>(_sensor.basis_uvw) * var_sensor_origin.x() + get<1>(_sensor.basis_uvw) * var_sensor_origin.y() + get<2>(_sensor.basis_uvw) * var_sensor_origin.z(),
                    get<0>(_sensor.basis_uvw) * var_direction.x() + get<1>(_sensor.basis_uvw) * var_direction.y() + get<2>(_sensor.basis_uvw) * var_direction.z(), vector3d<decimal>(var_factor_direction_z * _sensor.sensitivity), 1
                };
            }
            
            else if (param_lnm_type & lnm::non_specular_material)
            {
                const auto* var_material_map_kd = _material->map_kd();
                const vector3d<decimal> var_kd = (var_material_map_kd ? var_material_map_kd->eval(param_surface_point.tex_coord()) : _material->kd());
                decimal var_kd_max = var_kd.coord_max();
                decimal var_ks_max = _material->ks().coord_max();
                if (var_kd_max == 0 && var_ks_max == 0)
                {
                    var_kd_max = 1;
                    var_ks_max = 0;
                }
                decimal var_kd_ks = var_kd_max + var_ks_max;
                var_kd_max = var_kd_max / var_kd_ks;
                var_ks_max = var_ks_max / var_kd_ks;
                param_lnm_type = param_rng.uniform_real() < var_kd_max ? lnm::diffuse_material : lnm::glossy_material;
                
                // diffuse material direction sample evaluation
                if (param_lnm_type & lnm::diffuse_material)
                {
                    auto [var_normal, var_tangent, var_bi_tangent] = param_surface_point.ortho_normal_basis_wi(param_wi);
                    bool var_use_alpha = var_material_map_kd && !var_material_map_kd->alphas().empty();
                    decimal var_alpha_prob = var_use_alpha ? var_material_map_kd->eval(param_surface_point.tex_coord(), 1).x() : 1;
                    vector3d<decimal> var_cosine_weighted_direction = param_rng.cosine_weighted_direction();
                    direction_sample var_direction_sample_transparent_mask = {lnm::transparent_mask, param_surface_point.position(), -param_wi, vector3d<decimal>(1), var_kd_max};
                    direction_sample var_direction_sample_diffuse = direction_sample{lnm::diffuse_material, param_surface_point.position(), var_tangent * var_cosine_weighted_direction.x() + var_bi_tangent * var_cosine_weighted_direction.y() + var_normal * var_cosine_weighted_direction.z(), var_kd, var_kd_max};
                    return param_rng.uniform_real() > var_alpha_prob ? direction_sample{lnm::transparent_mask, param_surface_point.position(), -param_wi, vector3d<decimal>(1), var_kd_max} : direction_sample{lnm::diffuse_material, param_surface_point.position(), var_tangent * var_cosine_weighted_direction.x() + var_bi_tangent * var_cosine_weighted_direction.y() + var_normal * var_cosine_weighted_direction.z(), var_kd, var_kd_max};
                }
                
                // glossy material direction sample evaluation
                else if (param_lnm_type & lnm::glossy_material)
                {
                    auto [var_normal, var_tangent, var_bi_tangent] = param_surface_point.ortho_normal_basis_wi(param_wi);
                    const decimal var_theta = 2 * PI * param_rng.uniform_real();
                    const decimal var_phi = param_rng.uniform_real();
                    const vector3d<decimal> var_wh = normalize( (var_tangent*_material->roughness_x()*cos(var_theta)+var_bi_tangent*_material->roughness_y()*sin(var_theta)) * sqrt(var_phi/(1-var_phi)) + var_normal);
                    const vector3d<decimal> var_wo = reflect(param_wi, var_wh);
                    if (param_surface_point.op(param_wi, var_wo)){ return {}; }
                    return direction_sample{
                        lnm::glossy_material, param_surface_point.position(), var_wo, eval_bsdf(lnm::glossy_material, param_surface_point, param_wi, var_wo) / pdf_material(lnm::glossy_material, param_surface_point, param_wi, var_wo), var_ks_max};
                }
            }
            
            // mirror reflection
            else if (param_lnm_type & lnm::mirror_reflection)
            {
                return direction_sample{lnm::mirror_reflection, param_surface_point.position(), reflect(param_wi, param_surface_point.normal()), vector3d<decimal>(1), 1};
            }
            
            // fresnel reflection and refraction
            else if (param_lnm_type  & lnm::fresnel)
            {
                integer var_dot_wi_normal = dot(param_wi, param_surface_point.normal()) > 0;
                vector3d<decimal> var_normal = var_dot_wi_normal ? param_surface_point.normal() : -param_surface_point.normal();
                decimal var_eta = var_dot_wi_normal ? 1 / _material->ior() : _material->ior();
                auto var_wt = refract(param_wi, var_normal, var_eta);
                decimal var_fresnel = !var_wt ? 1 : [&]()
                {
                    decimal var_cosine = var_dot_wi_normal ? dot(param_wi, param_surface_point.normal()) : dot(*var_wt, param_surface_point.normal());
                    decimal var_ratio_eta = (1 - _material->ior()) / (1 + _material->ior());
                    var_ratio_eta = var_ratio_eta * var_ratio_eta;
                    return var_ratio_eta + (1 - var_ratio_eta) * pow(1 - var_cosine, 5);
                }();
                return param_rng.uniform_real() < var_fresnel ?
                    direction_sample{lnm::fresnel_reflection, param_surface_point.position(), reflect(param_wi, param_surface_point.normal()), vector3d<decimal>(1), 1} :
                direction_sample{lnm::fresnel_transmission, param_surface_point.position(), *var_wt, vector3d<decimal>(var_eta * var_eta), 1};
            }
            return {};
        }
        
        op<light_sample> sample_light(rng& param_rng,  geometry<>* param_geometry, const surface_point<> param_surface_point) const
        {
            if (_type & lnm::area_light)
            {
                const integer var_index_dist1d = _light_sampling_params.discrete_dist1d.sample(param_rng);
                const decimal var_sqrt_rng = sqrt(max(.0, param_rng.uniform_real()));
                const vector3d<decimal> var_vertex_a = param_geometry->positions[_face_indices[3 * var_index_dist1d].position];
                const vector3d<decimal> var_vertex_b = param_geometry->positions[_face_indices[3 * var_index_dist1d + 1].position];
                const vector3d<decimal> var_vertex_c = param_geometry->positions[_face_indices[3 * var_index_dist1d + 2].position];
                const surface_point<> var_surface_point_light(interpolate_barycentric<>(var_vertex_a, var_vertex_b, var_vertex_c, 1 - var_sqrt_rng, param_rng.uniform_real() * var_sqrt_rng), normalize(cross(var_vertex_b - var_vertex_a, var_vertex_c - var_vertex_a)), {});
                const vector3d<decimal> var_wo = var_surface_point_light.position() - param_surface_point.position();
                const vector3d<decimal> var_wo_normalized = normalize(var_wo);
                const decimal var_pdf_light = pdf_light(param_surface_point, var_surface_point_light, -var_wo_normalized);
                return var_pdf_light == 0 ? op<light_sample>{} : light_sample{var_wo_normalized, sqrt(dot(var_wo, var_wo)), eval_bsdf(lnm::area_light, var_surface_point_light, {}, -var_wo_normalized), var_pdf_light};
            }
            else if (_type & lnm::environment_light)
            {
                auto [var_u, var_v] = _env_light_sampling_params.discrete_dist2d.sample(param_rng);
                decimal var_theta = PI * var_v;
                decimal var_sin_theta = sin(var_theta);
                decimal var_phi = 2 * PI * var_u + _env_light_sampling_params.rotation;
                vector3d<decimal> var_wo(var_sin_theta * sin(var_phi), cos(var_theta), var_sin_theta * cos(var_phi));
                decimal var_pdf_light = pdf_light(param_surface_point, {}, -var_wo);
                return  var_pdf_light == 0 ? op<light_sample>{} : light_sample{var_wo, INFINITE, eval_bsdf(lnm::environment_light, {}, {}, -var_wo), var_pdf_light};
            }
            return {};
        }
    
        decimal pdf_light(const surface_point<>& param_surface_point, const surface_point<>& param_surface_point_light, const vector3d<decimal>& param_wo) const
        {
            if (_type & lnm::area_light)
            {
                decimal var_geometric_factor = geometric_factor(param_surface_point, param_surface_point_light);
                return var_geometric_factor == 0 ? 0 : _light_sampling_params.inverse_area / var_geometric_factor;
            }
            else if (_type & lnm::environment_light)
            {
                vector3d<> var_negative_wo = -param_wo;
                decimal var_arctan_wo = atan2(var_negative_wo.x(), var_negative_wo.z());
                var_arctan_wo = var_arctan_wo < 0 ? var_arctan_wo + 2 * PI : var_arctan_wo;
                decimal var_arctan_wo_rotated = (var_arctan_wo - _env_light_sampling_params.rotation) * .5 / PI;
                decimal var_u = var_arctan_wo_rotated - floor(var_arctan_wo_rotated);
                decimal var_v = acos(var_negative_wo.y()) / PI;
                decimal var_sqrt_wo_y = sqrt(1 - var_negative_wo.y() * var_negative_wo.y());
                return var_sqrt_wo_y == 0 ? 0 : _env_light_sampling_params.discrete_dist2d.pmf(var_u, var_v) / (2*PI*PI*var_sqrt_wo_y*abs(dot(-var_negative_wo, param_surface_point.normal())));
            }
            return 0;
        }
    
        vector3d<decimal> eval_bsdf(integer param_lnm_type, const surface_point<>& param_surface_point, const vector3d<>& param_wi, const vector3d<>& param_wo) const
        {
            if (param_lnm_type & lnm::area_light)
            {
                return dot(param_wo, param_surface_point.normal()) <= 0 ? vector3d<>() : _material->ke();
            }
            else if (param_lnm_type & lnm::environment_light)
            {
                const vector3d<> var_negative_wo = -param_wo;
                decimal var_arctan = atan2(var_negative_wo.x(), var_negative_wo.z());
                var_arctan = var_arctan < 0 ? var_arctan + 2 * PI : var_arctan;
                decimal var_tan = (var_arctan - _env_light_sampling_params.rotation) * .5 / PI;
                return _env_light_sampling_params.map.eval({var_tan - floor(var_tan), acos(var_negative_wo.y())/ PI, 0});
            }
            else if (param_lnm_type & lnm::glossy_material)
            {
                if (param_surface_point.op(param_wi, param_wo)) { return {}; }
                const auto* var_ptr_map_kd = _material->map_kd();
                vector3d<> var_diffuse = (var_ptr_map_kd ? var_ptr_map_kd->eval(param_surface_point.tex_coord()) : _material->kd()) * (1. / PI);
                const vector3d<> var_wh = normalize(param_wi + param_wo);
                const auto [var_normal, var_tangent, var_bi_tangent] = param_surface_point.ortho_normal_basis_wi(param_wi);
                const vector3d<> var_fresnel = _material->ks() + (vector3d<>(1) - _material->ks()) * pow(1 - dot(param_wo, var_wh), 5);
                return var_fresnel * (_material->aniso_ggx_distribution(var_wh, var_tangent, var_bi_tangent, var_normal) * _material->smith_g_aniso_ggx(param_wi, param_wo, var_tangent, var_bi_tangent, var_normal) / (4 * dot(param_wi, var_normal) * dot(param_wo, var_normal)));
            }
            else if (param_lnm_type & lnm::diffuse_material)
            {
                if (param_surface_point.op(param_wi, param_wo)) { return{}; }
                const auto* var_ptr_map_kd = _material->map_kd();
                const decimal var_alpha = (var_ptr_map_kd && !var_ptr_map_kd->alphas().empty()) ? var_ptr_map_kd->eval(param_surface_point.tex_coord(), 1).x() : 1;
                const auto var_diffuse = (var_ptr_map_kd ? var_ptr_map_kd->eval(param_surface_point.tex_coord()) : _material->kd()) * (var_alpha / PI);
                return var_diffuse;
            }
            return {};
        }
    
        decimal pdf_material(integer var_object_type, const surface_point<>& param_surface_point, const vector3d<decimal>& param_wi, const vector3d<decimal>& param_wo) const
        {
            return _material->pdf(var_object_type, param_surface_point, param_wi, param_wo);
        }
};

// obj file parsing utilities
bool check_space_character(CHAR param_char){ return param_char == ' ' || param_char == '\t'; };
bool check_command_token(CHAR*& param_token, const CHAR* param_command, integer param_index) 
{ return !strncmp(param_token, param_command, param_index) && check_space_character(param_token[param_index]); }
void ignore_spaces(CHAR *&param_token) { param_token += strspn(param_token, " \t"); }
void ignore_spaces_or_back_slash(CHAR *&param_token) { param_token += strcspn(param_token, "/ \t"); }

decimal parse_float(CHAR*& param_token)
{
    ignore_spaces(param_token);
    decimal var_value = atof(param_token);
    ignore_spaces_or_back_slash(param_token);
    return var_value;
}

vector3d<> parse_vector3d(CHAR*& param_token)
{
    vector3d<> var_vector3d;
    var_vector3d.set_x(parse_float(param_token));
    var_vector3d.set_y(parse_float(param_token));
    var_vector3d.set_z(parse_float(param_token));
    return var_vector3d;
}

integer parse_index(integer param_index, integer param_vertex_num)
{ 
    integer var_index = param_index < 0 ? param_vertex_num + param_index : param_index > 0 ? param_index - 1 : -1;
    return var_index; 
}

face_index parse_face_index(CHAR*& param_token, const geometry<>& param_geometry)
{
    face_index var_face_index;
    ignore_spaces(param_token);
    
    var_face_index.position = parse_index(atoi(param_token), integer(param_geometry.positions.size()));
    ignore_spaces_or_back_slash(param_token);
    if (param_token++ [0] != '/') { return var_face_index; }
    
    var_face_index.tex_coord = parse_index(atoi(param_token), integer(param_geometry.tex_coords.size()));
    ignore_spaces_or_back_slash(param_token);
    if (param_token++ [0] != '/') { return var_face_index; }
    
    var_face_index.normal = parse_index(atoi(param_token), integer(param_geometry.normals.size()));
    ignore_spaces_or_back_slash(param_token);
    return var_face_index;
}

void parse_string(CHAR*& param_token, CHAR param_name[]){ sscanf(param_token, "%s", param_name); }

class stage
{
    
    private:
        geometry<> _geometry;
        object _sensor;
        vector<object> _objects;
        vector<unique_ptr<texture>> _textures;
        vector<material<>> _materials;
        vector<integer> _light_indices;
        unordered_map<string, bool> _map_materials_lib;
        unordered_map<string, integer> _map_materials;
        unordered_map<string, integer> _map_textures;
        bool _env_light_enabled = false;
    
        // bvh data
        vector<bvh_node> _bvh_nodes;
        vector<triangle<>> _triangles;
        vector<integer> _triangle_indices;
    
    public:
        
        geometry<> get_geometry() { return _geometry;}
        geometry<>* get_geometry_ptr() { return &_geometry;}
        object* env_light() { return _env_light_enabled ? &_objects.front() : nullptr; }
        object* sensor_object_ptr() { return &_sensor;}
        object sensor_object() { return _sensor;}
        vector<integer> light_indices() const { return _light_indices;}
        tuple<const object&, decimal> sample_light(rng& param_rng) const
        {
            const integer var_total_light_indices = integer(_light_indices.size());
            const integer var_index = clamp(integer(param_rng.uniform_real() * var_total_light_indices), 0, var_total_light_indices - 1);
            return {_objects.at(_light_indices[var_index]), 1. / var_total_light_indices};
        }
    
        decimal pdf_light() { return 1. / _light_indices.size(); }
        void load_materials(string param_material_file_name)
        {
            ifstream var_file_material(param_material_file_name);
            const integer var_buffer_size = 4096;
            CHAR var_get_buffer_line[var_buffer_size], var_name[256];
            while (var_file_material.getline(var_get_buffer_line, var_buffer_size))
            {
                auto* var_token = var_get_buffer_line;
                ignore_spaces(var_token);
                if (check_command_token(var_token, "newmtl", 6))
                {
                    parse_string(var_token += 7, var_name);
                    _map_materials[var_name] = integer(_materials.size());
                    _materials.emplace_back();
                    continue;
                }
                if (_materials.empty()){ continue; }
                material<>& var_material = _materials.back();
                if (check_command_token(var_token, "Kd", 2)) { var_material.kd(parse_vector3d(var_token += 3)); }
                else if (check_command_token(var_token, "Ks", 2)) { var_material.ks(parse_vector3d(var_token += 3)); }
                else if (check_command_token(var_token, "Ke", 2))
                {
                    var_material.ke(parse_vector3d(var_token += 3));
                    var_material.type_ |= var_material.ke().coord_max() > 0 ? lnm::lights : 0;
                    var_material.type(var_material.type_);
                }
                else if (check_command_token(var_token, "illum", 5))
                {
                    ignore_spaces(var_token += 6);
                    const integer var_type = atoi(var_token);
                    var_material.type_ |= var_type == 7 ? lnm::fresnel : var_type == 5 ? lnm::mirror_reflection : lnm::non_specular_material;
                    var_material.type(var_material.type_);
                }
                else if (check_command_token(var_token, "map_Kd", 6))
                {
                    parse_string(var_token += 7, var_name);
                    auto iterator_map_textures = _map_textures.find(var_name);
                    if (iterator_map_textures != _map_textures.end())
                    {
                        var_material.map_kd(_textures[iterator_map_textures->second].get());
                        continue;
                    }
                    _map_textures[var_name] = integer(_textures.size());
                    _textures.emplace_back(new texture());
                    string var_file_texture_name = (file_system::path(param_material_file_name).remove_filename() / var_name).string();
                    _textures.back()->name(var_file_texture_name);
                    _textures.back()->loadppm(path(var_file_texture_name));
                    var_material.map_kd(_textures.back().get());
                    integer var_tex_nums = _textures.size();
                }
                else if (check_command_token(var_token, "Ni", 2)) { var_material.ior(parse_float(var_token += 3)); }
                else if (check_command_token(var_token, "Ns", 2)) { var_material.specular_exponent(parse_float(var_token += 3)); }
                else if (check_command_token(var_token, "aniso", 5)) { var_material.anisotropy(parse_float(var_token += 5)); }
                else { continue; }
            }
        }
    
        void load_obj(string param_obj_file_name, string param_env_map_file_name, decimal param_env_map_rotation)
        {
            cout << "\nstarted loading file: " << param_obj_file_name << "\n";
            const integer var_buffer_size = 4096;
            CHAR var_get_buffer_line[var_buffer_size], var_name[256];
            ifstream var_file_obj(param_obj_file_name);
            if (!param_env_map_file_name.empty())
            {
                _env_light_enabled = true;
                _objects.push_back({lnm::environment_light});
                _objects.back().init_env_light(param_env_map_file_name, param_env_map_rotation);
                _light_indices.push_back(0);
            }
            
            while (var_file_obj.getline(var_get_buffer_line, var_buffer_size))
            {
                CHAR* var_token = var_get_buffer_line;
                ignore_spaces(var_token);
                const integer var_total_objects = integer(_objects.size());
                if (check_command_token(var_token, "v", 1)) { _geometry.positions.emplace_back((parse_vector3d(var_token += 2))); }
                else if (check_command_token(var_token, "vn", 2)) { _geometry.normals.emplace_back((parse_vector3d(var_token += 3))); }
                else if (check_command_token(var_token, "vt", 2)) { _geometry.tex_coords.emplace_back((parse_vector3d(var_token += 3))); }
                else if (check_command_token(var_token, "f", 1))
                {
                    var_token += 2;
                    if (_materials.empty())
                    {
                        integer var_type = lnm::diffuse_material;
                        _materials.push_back({lnm::diffuse_material, 1});
                        _objects.push_back({lnm::diffuse_material, &_materials.back()});
                    }
                    
                    face_index var_face_index[4];
                    for (auto& index : var_face_index) 
                    { index = parse_face_index(var_token, _geometry); }
                    auto& var_ref_geometry_face_indices = _objects.back().face_indices();
                    var_ref_geometry_face_indices.insert(var_ref_geometry_face_indices.end(), {var_face_index[0], var_face_index[1], var_face_index[2]});
                    if (var_face_index[3].position != -1)
                    {
                        var_ref_geometry_face_indices.insert(var_ref_geometry_face_indices.end(), {var_face_index[0], var_face_index[2], var_face_index[3]});
                    }
                }
                else if (check_command_token(var_token, "usemtl", 6))
                {
                    var_token += 7;
                    parse_string(var_token, var_name);
                    auto& var_material = _materials[_map_materials.at(var_name)];

                    if(!var_material.type())
                    {
                        var_material.type(lnm::non_specular_material);
                    }
                    _objects.push_back({var_material.type(), &var_material});

                    if (var_material.type() & lnm::lights)
                    {
                         _light_indices.push_back(var_total_objects); 
                    }
                }
                else if (check_command_token(var_token, "mtllib", 6))
                {
                    parse_string(var_token += 7, var_name);
                    if(!_map_materials_lib[var_name])
                    {
                        load_materials(path((file_system::path(param_obj_file_name).remove_filename() / var_name).string()));
                        _map_materials_lib[var_name] = true;
                    }
                }
                else
                    continue;
            }
            
            for ( auto& var_object : _objects)
            {
                if (var_object.get_material() && var_object.get_material()->type() & lnm::area_light)
                {
                    var_object.init_area_light(_geometry);
                }
            }
            
            for ( auto& var_material : _materials)
            {
                if ( !(var_material.type() & lnm::non_specular_material)) { continue; }
                decimal var_radius = 2 / (2 + var_material.specular_exponent());
                decimal var_sqrt_anisotropy = sqrt(1 - var_material.anisotropy() * .9);
                var_material.roughness_x(max(1e-3, var_radius / var_sqrt_anisotropy));
                var_material.roughness_y(max(1e-3, var_radius / var_sqrt_anisotropy));
            }            
            cout << "\nfinished loading" << "\n";
            // print stats
            // cout << "\n total objects = " << _objects.size()<< "\n";
            // cout << "\n total materials = " << _materials.size()<< "\n";
            // cout << "\n total geometry positions = " << _geometry.positions.size() << "\n";
            // cout << "\n total geometry normals = " << _geometry.normals.size() << "\n";
            // cout << "\n total geometry tex coords = " << _geometry.tex_coords.size() << "\n";
        }
    
        void build_bvh()
        {
            cout << "\n building bvh" << "\n";

            // gather the triangles
            for (size_t var_object_index = 0; var_object_index < _objects.size(); var_object_index++)
            {
                auto& var_face_indices = _objects[var_object_index].face_indices();
                integer var_id_triangle = 0;
                for (size_t face_index = 0; face_index < var_face_indices.size(); face_index += 3)
                {
                    const vector3d<> var_a = _geometry.positions[var_face_indices[face_index].position];
                    const vector3d<> var_b = _geometry.positions[var_face_indices[face_index + 1].position];
                    const vector3d<> var_c = _geometry.positions[var_face_indices[face_index + 2].position];
                    _triangles.emplace_back(var_id_triangle, var_a, var_b, var_c, integer(var_object_index), integer(face_index));
                    integer current_size = _triangles.size();
                    var_id_triangle++;
                }
            }            
            const integer var_total_triangles = integer(_triangles.size());
            queue<tuple<integer, integer, integer>> traversal_queue;    //<node_index, start, end>
            traversal_queue.push({0, 0, var_total_triangles});              // start queue with root node
            _bvh_nodes.assign(2 * var_total_triangles -1, {});
            _triangle_indices.assign(var_total_triangles, 0);
            iota(_triangle_indices.begin(), _triangle_indices.end(), 0);
            mutex var_mutex;                                            // parallel queue
            condition_variable var_condition_variable;
            atomic<integer> var_triangles_processed = 0;
            atomic<integer> var_total_current_nodes = 1;
            bool finished = 0;
            
            auto task_bvh_create = [&]()
            {
                while (!finished)
                {
                    // create nodes for triangles within start and end
                    auto [var_node_index, var_start, var_end] = [&]() -> tuple<integer, integer, integer>
                    {
                        unique_lock<mutex> var_unique_lock_mutex(var_mutex);
                        if (!finished && traversal_queue.empty())
                        {
                            var_condition_variable.wait(var_unique_lock_mutex, [&]() { return finished || !traversal_queue.empty(); });
                        }
                        
                        if (finished) { return {}; }
                        auto var_queue = traversal_queue.front();
                        traversal_queue.pop();
                        return var_queue;
                    }();
                    
                    if (finished){ break; }
                    
                    // create bound for the node
                    bvh_node& var_bvh_node = _bvh_nodes[var_node_index];
                    if (var_node_index > _bvh_nodes.size())
                    {
                        cout << "\n out of bound access for index = " << var_node_index;
                    }
                    for (integer var_index = var_start; var_index < var_end; var_index++)
                    {
                        var_bvh_node.bound = merge_aabb(var_bvh_node.bound, _triangles[_triangle_indices[var_index]].bound());
                    }
                    
                    // sort triangles along the given axis
                    auto sort_triangles = [&, var_start = var_start, var_end = var_end](integer var_axis)
                    {
                        auto var_compare = [&](integer var_index1, integer var_index2)
                        {
                            return _triangles[var_index1].bound_centroid()[var_axis] < _triangles[var_index2].bound_centroid()[var_axis];
                        };
                        sort(&_triangle_indices[var_start], &_triangle_indices[var_end - 1] + 1, var_compare);
                    };
                    
                    // create leaf node
                    auto leaf_node = [&, var_start = var_start, var_end = var_end]()
                    {
                        var_bvh_node.leaf = true;
                        var_bvh_node.triangle_indices_start = var_start;
                        var_bvh_node.triangle_indices_end = var_end;
                        var_triangles_processed += var_end - var_start;
                        
                        if (var_triangles_processed == integer(_triangles.size()))
                        {
                            unique_lock<mutex> var_unique_lock_mutex(var_mutex);
                            finished = 1;
                            var_condition_variable.notify_all();
                        }
                    };
                    
                    // if net triangles is 1 then add leaf node
                    if (var_end - var_start < 2)
                    {
                        leaf_node();
                        continue;
                    }
                    
                    // SAH based split axis selection
                    decimal var_boundary = INFINITE;
                    integer var_split_index, var_split_axis;
                    for (integer var_axis = 0; var_axis < 3; var_axis++)
                    {
                        thread_local vector<decimal> var_total_triangles_left(var_total_triangles + 1);
                        thread_local vector<decimal> var_total_triangles_right(var_total_triangles + 1);
                        sort_triangles(var_axis);
                        aabb<> var_bound_left, var_bound_right;
                        
                        for (integer var_index = 0; var_index <= var_end - var_start; var_index++)
                        {
                            integer var_remaining = var_end - var_start - var_index;
                            var_total_triangles_left[var_index] = var_bound_left.surface_area() * var_index;
                            var_total_triangles_right[var_remaining] = var_bound_right.surface_area() * var_index;
                            var_bound_left = var_index < var_end - var_start ? merge_aabb(var_bound_left, _triangles[_triangle_indices[var_start + var_index]].bound()) : var_bound_left;
                            var_bound_right = var_remaining > 0 ? merge_aabb(var_bound_right, _triangles[_triangle_indices[var_start + var_remaining - 1]].bound()) : var_bound_right;
                        }
                        
                        for (integer var_index = 0; var_index < var_end - var_start; var_index++)
                        {
                            decimal var_sah_cost = 1 + (var_total_triangles_left[var_index] + var_total_triangles_right[var_index]) / var_bvh_node.bound.surface_area();
                            if (var_sah_cost < var_boundary)
                            {
                                var_boundary = var_sah_cost;
                                var_split_index = var_index;
                                var_split_axis  = var_axis;
                            }
                        }
                    }
                    if (var_boundary > var_end - var_start)
                    {
                        leaf_node();
                        continue;
                    }
                    
                    sort_triangles(var_split_axis);
                    integer var_between_index = var_start + var_split_index;
                    unique_lock<mutex> var_unique_lock_mutex(var_mutex);
                    traversal_queue.push({var_bvh_node.index_child_left = var_total_current_nodes++, var_start, var_between_index});
                    traversal_queue.push({var_bvh_node.index_child_right = var_total_current_nodes++, var_between_index, var_end});
                    var_condition_variable.notify_one();
                }
            };
            vector<thread> var_threads(omp_get_max_threads());
            for (auto& var_thread : var_threads){ var_thread = thread(task_bvh_create); }
            for (auto& var_thread : var_threads){ var_thread.join(); }            
            
            cout << "\n finished building " << "\n";
            cout << "\n total bvh nodes built = " << _bvh_nodes.size() << "\n";
        }
    
        op<hit_point> intersect_ray(const ray<>& param_ray, decimal param_hit_distance_low = EPSILON, decimal param_hit_distance_high = INFINITE)
        {
            op<hit_point_triangle> var_hit_point_triangle, var_hit_point_triangle_current;
            integer var_hit_index, var_total[99]{}, var_index = 0;
            while (var_index >= 0)
            {
                auto& var_bvh_node = _bvh_nodes.at(var_total[var_index--]);
                if (!var_bvh_node.bound.intersect_ray(param_ray, param_hit_distance_low, param_hit_distance_high)){ continue; }
                if (!var_bvh_node.leaf)
                {
                    var_total[++var_index] = var_bvh_node.index_child_left;
                    var_total[++var_index] = var_bvh_node.index_child_right;
                    continue;
                }
                
                for (integer var_index1 = var_bvh_node.triangle_indices_start; var_index1 < var_bvh_node.triangle_indices_end; var_index1++)
                {
                    if (var_hit_point_triangle_current = _triangles[_triangle_indices[var_index1]].intersect_ray(param_ray, param_hit_distance_low, param_hit_distance_high))
                    {
                        var_hit_point_triangle = var_hit_point_triangle_current;
                        param_hit_distance_high = var_hit_point_triangle_current->distance;
                        var_hit_index = var_index1;
                    }
                }
            }
            if (!var_hit_point_triangle)
            {
                return {};                
            }
            triangle<>& var_triangle = _triangles[_triangle_indices[var_hit_index]];
            object& var_object = _objects[var_triangle.object_index()];
            const vector3d<> var_point = param_ray.origin + param_ray.direction * param_hit_distance_high;
            const decimal var_u = var_hit_point_triangle->barycentric_u, var_v = var_hit_point_triangle->barycentric_v;
            const face_index& var_face_index1 = var_object.face_indices()[var_triangle.face_index()];
            const face_index& var_face_index2 = var_object.face_indices()[var_triangle.face_index() + 1];
            const face_index& var_face_index3 = var_object.face_indices()[var_triangle.face_index() + 2];
            const vector3d<> var_normal = var_face_index1.normal < 0  ? var_triangle.normal() : normalize(interpolate_barycentric(_geometry.normals[var_face_index1.normal], _geometry.normals[var_face_index2.normal], _geometry.normals[var_face_index3.normal], var_u, var_v));
            const vector3d<> var_tex_coord = var_face_index1.tex_coord < 0 ? 0 : interpolate_barycentric(_geometry.tex_coords[var_face_index1.tex_coord], _geometry.tex_coords[var_face_index2.tex_coord], _geometry.tex_coords[var_face_index3.tex_coord], var_u, var_v);
            return hit_point{{var_point, var_normal, var_tex_coord}, &var_object};
        }
};

//decimal clamp(double param_value){ return param_value<0 ? 0 : param_value>1 ? 1 : param_value; }
//decimal clamp(double param_value, double param_value_low, double param_value_high){ return param_value<param_value_low ? param_value_low : param_value>param_value_high ? param_value_high : param_value; }    
int toInt(double param_value){ return int(pow(clamp(param_value, 0.0, 1.0),1/2.2)*255+.5); }

//compile command <g++-9 -std=c++17 -o platinum_pt -lstdc++fs -fopenmp platinum_pt.cpp>
//run command <./platinum_pt.out ./salle_de_bain/render.args>
int main(integer param, CHAR **param_args)
{
    // parse render.args
    string var_file_render_args = param_args[1];
    if(!file_system::exists(var_file_render_args))
    {
        cout << "\n" << "no render args file supplied: aborting\n";
        return 0;
    }
    ifstream var_stream_render_ini;
    var_stream_render_ini.open(var_file_render_args.c_str());
    string var_key, var_value;
    map<string , string> var_map_render_args;
    while(var_stream_render_ini >> var_key >> var_value)
    {
        var_map_render_args[var_key] = var_value;
    }

    // setup render and camera parameters
    const auto obj_file_name = var_map_render_args["render_mesh"], env_map_file_name = var_map_render_args["env_map_file_name"];
    if(!file_system::exists(obj_file_name))
    {
        cout << "\n" << "no render mesh found: aborting\n";
        return 0;
    }
    cout << "\nfile to render : " << obj_file_name ;
    const auto lens_file_name = var_map_render_args["lens_file"];
    const auto samples_per_pixel = atoi(var_map_render_args["samples_per_pixel"].c_str()), max_path_length = atoi(var_map_render_args["max_path_length"].c_str());
    const auto russian_roulette_length = atoi(var_map_render_args["russian_roulette_length"].c_str());
    const auto var_radiance_clamp_max = atof(var_map_render_args["radiance_clamp_max"].c_str());
    const auto rotation_env_map = atoi(var_map_render_args["rotation_env_map"].c_str());
    const auto output_render_width = atoi(var_map_render_args["output_render_width"].c_str());
    const auto output_render_height = atoi(var_map_render_args["output_render_height"].c_str());
    const auto camera_position = vector3d<>(atof(var_map_render_args["camera_position_x"].c_str()), atof(var_map_render_args["camera_position_y"].c_str()), atof(var_map_render_args["camera_position_z"].c_str()));
    const auto camera_look_at = vector3d<>(atof(var_map_render_args["camera_look_at_x"].c_str()), atof(var_map_render_args["camera_look_at_y"].c_str()), atof(var_map_render_args["camera_look_at_z"].c_str()));
    const auto camera_up_vector = vector3d<>(atof(var_map_render_args["camera_up_vector_x"].c_str()), atof(var_map_render_args["camera_up_vector_y"].c_str()), atof(var_map_render_args["camera_up_vector_z"].c_str()));
    const auto fov_vertical = atof(var_map_render_args["fov_vertical"].c_str());
    const auto focus_distance = atof(var_map_render_args["focus_distance"].c_str());
    const auto sensor_diagonal_length = atof(var_map_render_args["sensor_diagonal_length"].c_str());
    const auto sensor_sensitivity = atof(var_map_render_args["sensor_sensitivity"].c_str());
    const auto aspect_ratio = decimal(output_render_width) / output_render_height;
    const auto output_driver = var_map_render_args["output_driver"];
    const auto output_file = var_map_render_args["output_file"];
    vector<vector3d<>> rgb_buffer_image(output_render_width * output_render_height, vector3d<>(0));
    
    // load obj file,environment map file(optional), build bvh and initialize sensor
    stage var_stage;
    var_stage.load_obj(obj_file_name, env_map_file_name, rotation_env_map);
    var_stage.build_bvh();    
    var_stage.sensor_object_ptr()->init_sensor(lens_file_name, camera_position, camera_look_at, camera_up_vector, fov_vertical, focus_distance, sensor_diagonal_length, sensor_sensitivity, aspect_ratio);
    
    // kajiya rendering integral equation
    auto lo = [&](rng& param_rng, const vector3d<>& pixel_position)
    {
        vector3d<> var_radiance_output, var_throughput(1), var_wi = pixel_position;
        integer var_object_type_current = lnm::sensor;
        object* var_object_current = var_stage.sensor_object_ptr();
        surface_point<> var_surface_point_current;
        
        for (integer var_path_length = 0; var_path_length < max_path_length; var_path_length++)
        {
            integer obj_type = lnm::non_specular_material;
            const bool next_event_estimated = !var_stage.light_indices().empty() && var_path_length > 0 && var_object_type_current & lnm::non_specular_material;            
            auto var_direction_sample = var_object_current->sample_direction(param_rng, var_object_type_current, var_surface_point_current, var_wi  );

            if (!var_direction_sample)
            { 
                break; 
            }
            auto [var_object_type, var_ray_sampled, var_evaluated_w, var_prob_component_select] = *var_direction_sample;
            var_throughput = var_throughput / var_prob_component_select;
            
            if (next_event_estimated)
            {
                auto [var_object_light, var_inverse_light_indices] = var_stage.sample_light(param_rng);
                auto var_light_sample = var_object_light.sample_light(param_rng, var_stage.get_geometry_ptr(), var_surface_point_current);
                if (var_light_sample)
                {
                    auto [var_wo, var_sample_distance, var_evaluated_le, var_evaluated_prob] = *var_light_sample;
                    if (!var_stage.intersect_ray({var_surface_point_current.position(), var_wo}, EPSILON, var_sample_distance * (1 - EPSILON)))
                    {
                        const auto var_eval_bsdf = var_object_current->eval_bsdf(var_object_type, var_surface_point_current, var_wi, var_wo) * var_evaluated_le;
                        const auto var_eval_pdf = (var_object_current->pdf_material(var_object_type, var_surface_point_current, var_wi, var_wo) + var_evaluated_prob * var_inverse_light_indices);
                        var_radiance_output = var_radiance_output + var_throughput * var_eval_bsdf / var_eval_pdf;
                    }
                }
            }
            
            // update throughput
            var_throughput = var_throughput * var_evaluated_w;
            
            // continue tracing
            auto var_hit_point = var_stage.intersect_ray(var_ray_sampled);
            if (!var_hit_point)
            {
                // we might hit environment light otherwise break
                var_hit_point = hit_point{{}, var_stage.env_light()};
                if (!var_hit_point->object_hit_ptr)
                {
                    break;                    
                }
            }
            
            // accumulate light contribution
            auto [var_surface_point_hit, var_object_hit_ptr] = *var_hit_point;
            integer var_object_hit_type = var_object_hit_ptr->type();
            if (var_object_hit_type & lnm::lights)
            {
                const auto var_eval_bsdf = var_object_hit_ptr->eval_bsdf(var_object_hit_type & lnm::lights, var_surface_point_hit, {}, -var_ray_sampled.direction);
                const auto var_pdf_light = var_object_hit_ptr->pdf_light(var_surface_point_current, var_surface_point_hit, -var_ray_sampled.direction);
                const auto var_pdf_material = var_object_current->pdf_material(var_object_hit_type, var_surface_point_current, var_wi, var_ray_sampled.direction);
                const auto var_eval_pdf = (var_path_length == 0 || var_object_type & lnm::specular_material ? 1 :  var_pdf_light * var_stage.pdf_light() / var_pdf_material + 1);
                var_radiance_output = var_radiance_output + var_throughput * var_eval_bsdf / var_eval_pdf;
            }
            
            // russian roulette termination
            if (var_path_length > russian_roulette_length)
            {
                decimal var_termination_prob_threshold = max(.2, 1 - var_throughput.coord_max());
                // boost throughput by continuation pdf
                var_throughput = var_throughput / (1 - var_termination_prob_threshold);
                if (param_rng.uniform_real() < var_termination_prob_threshold) 
                { 
                    break; 
                }
            }
            
            // update path ray intersection info
            var_wi = -var_ray_sampled.direction;
            var_object_current = var_object_hit_ptr;
            var_surface_point_current = var_surface_point_hit;
            var_object_type_current = var_object_current->type() & ~lnm::lights;
        }        
        return var_radiance_output;
    };
    
    cout << "\nrendering started.....\n";
    cout << "\nimage width is = " << output_render_width;
    cout << "\nimage height is = " << output_render_height;
    cout << "\nsamples per pixel is = " << samples_per_pixel << "\n";
    
    //integer h = output_render_height, w = output_render_width;
    #pragma omp parallel for schedule(dynamic, 1)
    for (integer index_pixel = 0; index_pixel < output_render_width * output_render_height; index_pixel++)
    {
        thread_local rng var_rng(137 + omp_get_thread_num());
        for (integer index_sample = 0; index_sample < samples_per_pixel; index_sample++)
        {
            const vector3d<> pixel_position((index_pixel%output_render_width+var_rng.uniform_real())/output_render_width, (index_pixel/output_render_width+var_rng.uniform_real())/output_render_width, 0);
            vector3d<> var_rad_out = lo(var_rng, pixel_position);
            var_rad_out = vector3d<>(clamp(var_rad_out.x(),0.0, var_radiance_clamp_max), clamp(var_rad_out.y(), 0.0, var_radiance_clamp_max), clamp(var_rad_out.z(), 0.0, var_radiance_clamp_max));
            rgb_buffer_image[index_pixel] = rgb_buffer_image[index_pixel] + var_rad_out / samples_per_pixel;
        }
    }
    cout << "\nrendering over";
    
    // save buffer image to file    
    if(output_driver.compare("ppm") == 0)
    {       
        string var_file_name = output_file + "." + output_driver; 
        FILE *var_file = fopen(var_file_name.c_str(), "w");
        fprintf(var_file, "P3\n%d %d\n%d\n", output_render_width, output_render_height, 255); 
        for (int index_pixel=0; index_pixel<output_render_width*output_render_height; index_pixel++)
        { 
            fprintf(var_file,"%d %d %d ", toInt(rgb_buffer_image[index_pixel].x()), toInt(rgb_buffer_image[index_pixel].y()), toInt(rgb_buffer_image[index_pixel].z()));
        }            
        fclose(var_file);
        cout << "\n" << "image file saved : " <<  var_file_name << "\n\n";
    }
    else
    {
        cout << "\n" << "unsupported driver requested" << "\n" << "Only PPM files are supported";
    }   

    return 0;
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
