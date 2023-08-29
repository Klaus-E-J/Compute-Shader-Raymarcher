#ifndef COMPUTE_SHADER_H
#define COMPUTE_SHADER_H

#include <GL/glew.h>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class ComputeShader
{
public:
    GLuint computeProgram;
    
    ComputeShader(const char* computePath)
    {
        std::string computeCode;
        
        std::ifstream cShaderFile;
        // ensure ifstream objects can throw exceptions:
        cShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
        
        try
        {
            // Open file
            cShaderFile.open(computePath);
            std::stringstream cShaderStream;
            // read file's buffer contents into streams
            cShaderStream << cShaderFile.rdbuf();
            // Close file handlers
            cShaderFile.close();
            computeCode = cShaderStream.str();
        }
        catch (std::ifstream::failure& e)
        {
            std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ: " << e.what() << std::endl;
        }
        
        const char* cShaderCode = computeCode.c_str();
        
        GLuint compute;
        
        // compute Shader
        compute = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(compute, 1, &cShaderCode, NULL);
        glCompileShader(compute);
        checkCompileErrors(compute, "COMPUTE");
        
        computeProgram = glCreateProgram();
        
        glAttachShader(computeProgram, compute);
        glLinkProgram(computeProgram);
        
        checkCompileErrors(computeProgram, "COMPUTE PROGRAM");
        
        glDeleteShader(compute);
    }
    
    void Use(int x_resolution, int y_resolution)
    {
        glUseProgram(computeProgram);
        glDispatchCompute(ceil(x_resolution / 8), ceil(y_resolution / 8), 1); // IMPORTANT: Allways use these two lines before using the conventional shader program
		glMemoryBarrier(GL_ALL_BARRIER_BITS);                // or else, the compute program will not work
    }
    
    void Delete()
    {
        glDeleteProgram(computeProgram);
    }
    
private:
    // utility function for checking shader compilation/linking errors.
    // ------------------------------------------------------------------------
    void checkCompileErrors(unsigned int shader, std::string type)
    {
        int success;
        char infoLog[1024];
        if (type != "PROGRAM")
        {
            glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
            if (!success)
            {
                glGetShaderInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::SHADER_COMPILATION_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
        else
        {
            glGetProgramiv(shader, GL_LINK_STATUS, &success);
            if (!success)
            {
                glGetProgramInfoLog(shader, 1024, NULL, infoLog);
                std::cout << "ERROR::PROGRAM_LINKING_ERROR of type: " << type << "\n" << infoLog << "\n -- --------------------------------------------------- -- " << std::endl;
            }
        }
    }
    
};
#endif
