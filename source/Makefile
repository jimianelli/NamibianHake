
ifdef ComSpec
    RM=del /F /Q
    COPY=copy
else
    RM=rm -rf
    COPY=cp
endif

all: nh 

nh: nh.tpl 
	@admb -f nh.tpl 


clean:
	@$(RM) nh.cpp
	@$(RM) nh.htp
	@$(RM) nh.obj
