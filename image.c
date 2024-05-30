#include "lodepng.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14
#define radius 1
#define sigma 0.4
#define epsilon 100

typedef struct  Node {
	unsigned char r, g, b, a;
	struct Node *up, *down, *left, *right, *parent;
	int rank;
}Node;

float gaussian(int x, int y) {
	return exp(-(x * x + y * y) / (2 * sigma * sigma)) / (2 * PI * sigma * sigma);
}

void filter_gauss(unsigned char* picture, unsigned char* gauss_picture, unsigned w, unsigned h) {
	int kernel_size = 2 * radius + 1, i, j, x, y, kx, ky, pic_index, gauss_index;
	float* kernel = malloc(kernel_size * kernel_size * sizeof(float));
	float kernel_sum = 0.0, kernel_value, r, g, b;
	for (i = -radius; i <= radius; i++)
		for (j = -radius; j <= radius; j++) {
			kernel[(i + radius) * kernel_size + (j + radius)] = gaussian(i, j, sigma);
			kernel_sum += kernel[(i + radius) * kernel_size + (j + radius)];
		}

	for (i = 0; i < kernel_size * kernel_size; i++)
		kernel[i] /= kernel_sum;

	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			r = 0.0, g = 0.0, b = 0.0;
			for (i = -radius; i <= radius; i++)
				for (j = -radius; j <= radius; j++) {
					kx = x + j;
					ky = y + i;
					kx = (kx < 0) ? 0 : kx;
					ky = (ky < 0) ? 0 : ky;
					ky = (ky >= h) ? (h - 1) : ky;
					kx = (kx >= w) ? (w - 1) : kx;
					pic_index = (ky * w + kx) * 4;
					kernel_value = kernel[(i + radius) * kernel_size + (j + radius)];
					r += picture[pic_index] * kernel_value;
					g += picture[pic_index + 1] * kernel_value;
					b += picture[pic_index + 2] * kernel_value;
				}

			int gauss_index = (y * w + x) * 4;
			gauss_picture[gauss_index] = (r > 255.0) ? (unsigned char)255.0 : (unsigned char)r;
			gauss_picture[gauss_index + 1] = (g > 255.0) ? (unsigned char)255.0 : (unsigned char)g;
			gauss_picture[gauss_index + 2] = (b > 255.0) ? (unsigned char)255.0 : (unsigned char)b;
		}
	}
	free(kernel);
}

void components(int w, int h, unsigned char* picture, Node* nodes) {  //filling our array of nodes. Npw every node's parent is the same node.
	int y, x;
	for (y = 0; y < h; y++)
		for (x = 0; x < w; x++) {
			Node* current_node = &nodes[y * w + x];
			current_node->r = picture[(y * w + x) * 4];
			current_node->g = picture[(y * w + x) * 4 + 1];
			current_node->b = picture[(y * w + x) * 4 + 2];
			current_node->a = picture[(y * w + x) * 4 + 3];
			current_node->left = x > 0 ? &nodes[y * w + x - 1] : NULL;
			current_node->right = x < w - 1 ? &nodes[y * w + x + 1] : NULL;
			current_node->up = y > 0 ? &nodes[(y - 1) * w + x] : NULL;
			current_node->down = y < h - 1 ? &nodes[(y + 1) * w + x] : NULL;
			current_node->parent = current_node;
			current_node->rank = 0;
		}
}


Node* find_parent(Node* x) {  //for our union–find data structure
	if (x->parent != x)
		x->parent = find_parent(x->parent);
	return x->parent;
}

void union_set(Node* x, Node* y) { //x and y are nodes that we want to union
	if (x->r < 20 && y->r < 20)
		return;
	Node* xset = find_parent(x);
	Node* yset = find_parent(y);

	double color_difference = sqrt(pow(x->r - y->r, 2) * 3);
	if (xset != yset && color_difference < epsilon) {
		if (xset->rank > yset->rank)
			yset->parent = xset;
		else {
			xset->parent = yset;
			if (xset->rank == yset->rank)
				yset->rank++;
		}
	}
}

void find_components(int w, int h, Node* nodes) {
	int x, y;
	for (y = 0; y < h; y++) {
		for (x = 0; x < w; x++) {
			Node* node = &nodes[y * w + x];
			if (node->left != NULL)
				union_set(node, node->left);
			if (node->right != NULL)
				union_set(node, node->right);
			if (node->up != NULL)
				union_set(node, node->up);
			if (node->down != NULL)
				union_set(node, node->down);
		}
	}
}


void colouring(int w, int h, Node* nodes, unsigned char* output_image) {
	int* component_sizes = calloc(w * h, sizeof(int));
	int total_components = 0, i;

	srand(time(NULL));
	for (i = 0; i < w * h; i++) {
		Node* p = find_parent(&nodes[i]);
		if (p == &nodes[i]) {
			if (component_sizes[i] < 15) {
				p->r = 0;
				p->g = 0;
				p->b = 0;
			}
			else {
				p->r = rand() % 256 + 100;
				p->g = rand() % 256 + 100;
				p->b = rand() % 256 + 100;
			}
			total_components++;
		}
		output_image[4 * i + 0] = p->r;
		output_image[4 * i + 1] = p->g;
		output_image[4 * i + 2] = p->b;
		output_image[4 * i + 3] = 255;
		component_sizes[p - nodes]++;
	}
	free(component_sizes);
}





char* load_png_file(const char* filename, int* width, int* height) {
	unsigned char* image = NULL;
	int error = lodepng_decode32_file(&image, width, height, filename);
	if (error) {
		printf("error %u: %s\n", error, lodepng_error_text(error));
		return NULL;
	}
	return (image);
}

void to_grayscale(char* picture, unsigned int width, unsigned int height) {
	char* p = picture;
	int y, x;
	unsigned char  maxR = 0, maxG = 0, maxB = 0; //for normalizing
	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			unsigned char  r = *p;
			unsigned char  g = *(p + 1);
			unsigned char  b = *(p + 2);
			unsigned char  a = *(p + 3);
			maxR = (r > maxR) ? r : maxR; //choosing maximum from all channels of all picture
			maxG = (g > maxG) ? g : maxG;
			maxB = (b > maxB) ? b : maxB;
			p += 4;  //moving to next pixel
		}

	p = picture;  //getting back

	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			unsigned char r = *p;
			unsigned char g = *(p + 1);
			unsigned char b = *(p + 2);
			unsigned char a = *(p + 3);
			unsigned char gray = (unsigned char)(255.0 * ((float)r / maxR + (float)g / maxG + (float)b / maxB) / 3); // average of r, g, and b
			*p = *(p + 1) = *(p + 2) = gray; //grayscale is when r, g, b are equal
			p += 4;
		}
	}
}

void sobel(unsigned char* picture, unsigned char* new_picture, unsigned width, unsigned height) {
	int dx, dy, grad, i, j;
	unsigned int x, y;

	int sobelX[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };
	int sobelY[3][3] = { {-1, -2, -1}, {0, 0, 0}, {1, 2, 1} };

	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++) {
			dx = 0;
			dy = 0;
			unsigned int rgb, index;
			if (y == 0 || y == height - 1 || x == 0 || x == width - 1)
				grad = 0; //for pixels on the boundary
			else {
				for (i = -1; i <= 1; i++)
					for (j = -1; j <= 1; j++) {
						index = ((x + i) + (y + j) * width) * 4;
						rgb = picture[index]; //we take red channel (all all channels are the same because black&white)
						dx += rgb * sobelX[i + 1][j + 1];
						dy += rgb * sobelY[i + 1][j + 1];
					}
				grad = sqrt(dx * dx + dy * dy); //calculating grad's magnitude
			}

			//normalizing
			grad = (grad > 255) ? 255 : grad;
			grad = (grad < 0) ? 0 : grad;
			index = 4 * (x + y * width);
			new_picture[index] = grad;  //all channels are the same because black&white
			new_picture[index + 1] = grad;
			new_picture[index + 2] = grad;
			new_picture[index + 3] = picture[index + 3]; //taking alpha channel from previous picture
		}
	return 0;

}

int main() {

	int w = 0, h = 0;
	int k = 0;

	char* filename = "input1.png";
	char* output_filename = "output.png";

	char* picture = load_png_file(filename, &w, &h);
	if (picture == NULL) {
		printf("I can't read the picture %s. Error.\n", filename);
		return -1;
	}


	to_grayscale(picture, w, h); //now picture is BW

	// memory for after-sobel image
	unsigned char* sobel_picture = malloc(w * h * 4);
	if (sobel_picture == NULL) {
		printf("Не удалось выделить память\n");
		free(picture);
		return 1;
	}

	//memory for gauss picture
	unsigned char* gauss_picture = malloc(w * h * 4);
	if (gauss_picture == NULL) {
		printf("Не удалось выделить память\n");
		free(picture);
		return 1;
	}

	//memory for nodes of the graph
	Node* nodes = malloc(sizeof(Node) * w * h);
	if (nodes == NULL) {
		printf("Не удалось выделить память\n");
		free(nodes);
		return 1;
	}

	//memory for final picture (after colouring)
	unsigned char* output_image = malloc(w * h * 4);
	if (output_image == NULL) {
		printf("Не удалось выделить память\n");
		free(output_image);
		return 1;
	}


	filter_gauss(picture, gauss_picture, w, h);
	sobel(gauss_picture, sobel_picture, w, h);
	components(w, h, sobel_picture, nodes);

	find_components(w, h, nodes);
	colouring(w, h, nodes, output_image);

	int error = lodepng_encode32_file(output_filename, output_image, w, h);  //saving new_picture
	if (error)
		printf("error %u: %s\n", error, lodepng_error_text(error));


	free(picture);
	free(gauss_picture);
	free(sobel_picture);
	free(nodes);
	return 0;
}
