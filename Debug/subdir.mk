################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ECEFNRT.c \
../rsr_utils.c 

O_SRCS += \
../ECEFNRT.o \
../rsr_utils.o 

OBJS += \
./ECEFNRT.o \
./rsr_utils.o 

C_DEPS += \
./ECEFNRT.d \
./rsr_utils.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


