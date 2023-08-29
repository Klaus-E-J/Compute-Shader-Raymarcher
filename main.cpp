#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <cmath>
#include <string>
#include <glm/glm.hpp>

#include "ShaderClass.h"
#include "ComputeShaderClass.h"

// Function definition that un run everytime the window resizes. This "normalizes" the aspect of it
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Function to process the window input
void processInput(GLFWwindow *window);


double prevTime = 0.0;
int frameCount = 0;
int fps = 0;

void calculateFPS(float currentTime) {
    frameCount++;

    if (currentTime - prevTime >= 1.0) {
        fps = frameCount;
        frameCount = 0;
        prevTime = currentTime;
    }
}

glm::dvec2 mouse_pos;

void LockMouse(GLFWwindow* window) {
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwSetCursorPos(window, windowWidth / 2, windowHeight / 2);
    mouse_pos.x = windowWidth / 2;
    mouse_pos.y = windowHeight / 2;
}

glm::dvec2 getMouseDelta(GLFWwindow* window) {
    // Get the current mouse position
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    glm::dvec2 mousePos(xpos, ypos);

    // Calculate the difference in mouse position since the last frame
    glm::dvec2 mouseDelta = mousePos - mouse_pos;

    // Update the locked mouse position to the center of the window
    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    glfwSetCursorPos(window, windowWidth / 2.0, windowHeight / 2.0);

    // Update the mouse position for the next frame
    mouse_pos.x = windowWidth / 2.0;
    mouse_pos.y = windowHeight / 2.0;

    // Return the mouse delta
    return mouseDelta;
}

// window settings
const unsigned int SCR_WIDTH = 1200;
const unsigned int SCR_HEIGHT = 1000;

const unsigned int LOW_SCR_WIDTH = 350;
const unsigned int LOW_SCR_HEIGHT = 350;

bool vSync = false;

GLfloat vertices[] =
{
	-1.0f, -1.0f , 0.0f, 0.0f, 0.0f,
	-1.0f,  1.0f , 0.0f, 0.0f, 1.0f,
	 1.0f,  1.0f , 0.0f, 1.0f, 1.0f,
	 1.0f, -1.0f , 0.0f, 1.0f, 0.0f,
};

GLuint indices[] =
{
	0, 2, 1,
	0, 3, 2
};

int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif


    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Compute", glfwGetPrimaryMonitor(), NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(vSync);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    //GLEW initialization
    glfwMakeContextCurrent(window);
    if ( glewInit() != GLEW_OK )
       std::cout << "error" << std::endl;
    
    // Do any OpenGL initialization, including VAO, VBO, EBO, Shader compilation, etc.
    
    GLuint VAO, VBO, EBO;
    glCreateVertexArrays(1, &VAO);
    glCreateBuffers(1, &VBO);
    glCreateBuffers(1, &EBO);
    
    glNamedBufferData(VBO, sizeof(vertices), vertices, GL_STATIC_DRAW);
    glNamedBufferData(EBO, sizeof(indices), indices, GL_STATIC_DRAW);
    
    glEnableVertexArrayAttrib(VAO, 0);
    glVertexArrayAttribBinding(VAO, 0, 0);
    glVertexArrayAttribFormat(VAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
    
    glEnableVertexArrayAttrib(VAO, 1);
	glVertexArrayAttribBinding(VAO, 1, 0);
	glVertexArrayAttribFormat(VAO, 1, 2, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat));

	glVertexArrayVertexBuffer(VAO, 0, VBO, 0, 5 * sizeof(GLfloat));
	glVertexArrayElementBuffer(VAO, EBO);


	GLuint screenTex;
	glCreateTextures(GL_TEXTURE_2D, 1, &screenTex);
	glTextureParameteri(screenTex, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTextureParameteri(screenTex, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTextureParameteri(screenTex, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTextureParameteri(screenTex, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTextureStorage2D(screenTex, 1, GL_RGBA32F, LOW_SCR_WIDTH, LOW_SCR_HEIGHT);
	glBindImageTexture(0, screenTex, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);

	Shader shaderSources("shaders/vertex.glsl", "shaders/fragment.glsl");
	ComputeShader computeSource("shaders/compute.glsl");
	
	glm::dvec3 pos(0.0, 0.0, 0.0);
	glm::dvec3 dir(0.0, 0.0, 0.0);
	constexpr double cstSpeed = 0.1;
	double moveSpeed = cstSpeed;
	
	// Create a placeholder variable to serve as a "link" betwen the uniform in the compute shader and the variable that we may use coming from the main 
	// application
	GLint uPosLocation = glGetUniformLocation(computeSource.computeProgram, "u_pos");
	GLint uResolutionLocation = glGetUniformLocation(computeSource.computeProgram, "u_resolution");
	GLint uDirLocation = glGetUniformLocation(computeSource.computeProgram, "u_dir");
	GLint uTimeLocation = glGetUniformLocation(computeSource.computeProgram, "u_time");
	
	/*
	GLuint screenVertexShader = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(screenVertexShader, 1, &screenVertexShaderSource, NULL);
	glCompileShader(screenVertexShader);
	
	GLuint screenFragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(screenFragmentShader, 1, &screenFragmentShaderSource, NULL);
	glCompileShader(screenFragmentShader);

	GLuint screenShaderProgram = glCreateProgram();
	glAttachShader(screenShaderProgram, screenVertexShader);
	glAttachShader(screenShaderProgram, screenFragmentShader);
	glLinkProgram(screenShaderProgram);

	glDeleteShader(screenVertexShader);
	glDeleteShader(screenFragmentShader);
    
    // Actual Compute shader
    GLuint computeShader = glCreateShader(GL_COMPUTE_SHADER);
	glShaderSource(computeShader, 1, &screenComputeShaderSource, NULL);
	glCompileShader(computeShader);

	GLuint computeProgram = glCreateProgram();
	glAttachShader(computeProgram, computeShader);
	glLinkProgram(computeProgram);
	*/
    
    // For finding work group limitations
    
    int work_grp_cnt[3];
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &work_grp_cnt[0]);
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &work_grp_cnt[1]);
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &work_grp_cnt[2]);
	std::cout << "Max work groups per compute shader" << 
		" x:" << work_grp_cnt[0] <<
		" y:" << work_grp_cnt[1] <<
		" z:" << work_grp_cnt[2] << "\n";

	int work_grp_size[3];
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &work_grp_size[0]);
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 1, &work_grp_size[1]);
	glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 2, &work_grp_size[2]);
	std::cout << "Max work group sizes" <<
		" x:" << work_grp_size[0] <<
		" y:" << work_grp_size[1] <<
		" z:" << work_grp_size[2] << "\n";

	int work_grp_inv;
	glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &work_grp_inv);
	std::cout << "Max invocations count per work group: " << work_grp_inv << "\n";
    
	LockMouse(window);
	while (!glfwWindowShouldClose(window))
    {
		float currentTime = glfwGetTime();
		
		//std::cout << "Time: " << currentTime << '\n';
		
		glm::ivec2 mouseDelta = getMouseDelta(window);
		//std::cout << "Mouse Delta X: " << mouseDelta.x << '\n';
		//std::cout << "Mouse Delta Y: " << mouseDelta.y << '\n';
		
		dir.x += float(mouseDelta.x) * 0.002;
		
		dir.y += float(mouseDelta.y) * 0.002;
		
		//Camera clamping
		if(dir.y <= -1.75)
		{
			dir.y = -1.75f;
		}
        
		if(dir.y >= 1.75)
		{
			dir.y = 1.75;
		}
		
		glm::dvec3 forward(std::sin(dir.x), 0.0, std::cos(dir.x));
		glm::dvec3 right(std::sin(dir.x - M_PI / 2), 0.0, std::cos(dir.x - M_PI / 2));
		
		
		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
		{
			pos += moveSpeed * forward;
		}
		if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
		{
			pos += moveSpeed * right;
		}
		if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
		{
			pos -= moveSpeed * forward;
		}
		if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
		{
			pos -= moveSpeed * right;
		}
		if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
		{
			pos.y += moveSpeed;
		}
		if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
		{
			pos.y -= moveSpeed;
		}
		if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
		{
			moveSpeed = 10.0f;
		}
		else if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_RELEASE)
		{
			moveSpeed = cstSpeed;
		}
		
		//std::cout << "X pos: " << pos.x << '\n';
        //std::cout << "Z pos: " << pos.z << '\n';
		
		//pos.y += 0.01;
		computeSource.Use(LOW_SCR_WIDTH, LOW_SCR_HEIGHT);
		// Set uniform values in the compute shader
		glUniform2f(uResolutionLocation, LOW_SCR_WIDTH, LOW_SCR_HEIGHT);
        glUniform3f(uPosLocation, pos.x, pos.y, pos.z);
		glUniform3f(uDirLocation, dir.x, dir.y, dir.z);
		glUniform1f(uTimeLocation, currentTime);
		shaderSources.Use();
		//glDispatchCompute(ceil(SCR_WIDTH / 8), ceil(SCR_HEIGHT / 8), 1);
		//glMemoryBarrier(GL_ALL_BARRIER_BITS);
		/*
        // To use the compute shader you need to put in inside the main loop
        glUseProgram(computeProgram);
		glDispatchCompute(ceil(SCR_WIDTH / 8), ceil(SCR_HEIGHT / 8), 1);
		glMemoryBarrier(GL_ALL_BARRIER_BITS);
        
        glUseProgram(screenShaderProgram);
        */
		glBindTextureUnit(0, screenTex);
		glUniform1i(glGetUniformLocation(shaderSources.ID, "screen"), 0);
		glBindVertexArray(VAO);
		glDrawElements(GL_TRIANGLES, sizeof(indices) / sizeof(indices[0]), GL_UNSIGNED_INT, 0);

        
        // input
        // -----
        processInput(window);
		
		calculateFPS(currentTime);
        std::string title = "FPS: " + std::to_string(fps);
        glfwSetWindowTitle(window, title.c_str());
        
        
        glfwSwapBuffers(window);
        glfwPollEvents();
		
    }
    computeSource.Delete();
    
}

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}
