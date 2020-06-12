# Build secure-gwas docker image
image_name = secure-gwas
build : 
	docker build -t $(image_name) .

# Run secure-gwas container
container_name = $(image_name)_$(shell date +'%s')
run : 
	docker run --name $(container_name) $(image_name)