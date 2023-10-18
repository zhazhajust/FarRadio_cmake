#include "FaradioMPI.hpp"

FaradioMPI::FaradioMPI() {
    this->init(0, NULL);
};

FaradioMPI::~FaradioMPI() {
    this->finalize();
};

std::shared_ptr<FaradioMPI> FaradioMPI::getFaradioMPI(){
    std::call_once(singletonFlag, [&](){
        //instance.reset(new FaradioMPI());
        instance = std::make_shared<FaradioMPI>();
    });
    return instance;
};

void FaradioMPI::init(int argc, char** argv) {
    int flag = 0;
    MPI_Initialized(&flag);
    if(!flag) MPI_Init(&argc, &argv);
    this->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(this->comm, (int*)&this->_rank);
    MPI_Comm_size(this->comm, (int*)&this->_size);
}

void FaradioMPI::finalize() {
    MPI_Finalize();
};

void FaradioMPI::send(double* data, int count, int dest, int tag) {
    MPI_Send(data, count, MPI_DOUBLE, dest, tag, this->comm);
};

void FaradioMPI::recv(double* data, int count, int source, int tag) {
    MPI_Recv(data, count, MPI_DOUBLE, source, tag, this->comm, MPI_STATUS_IGNORE);
};

void FaradioMPI::reduce(double* send_data, double* recv_data, int count, int root) {
    MPI_Reduce(send_data, recv_data, count, MPI_DOUBLE, MPI_SUM, root, this->comm);
};

void FaradioMPI::reduce(void* send_data, void* recv_data, int count, int root) {
    MPI_Reduce(send_data, recv_data, count, MPI_DOUBLE, MPI_SUM, root, this->comm);
};

void FaradioMPI::gather(double* send_data, double* recv_data, int count, int root) {
    MPI_Gather(send_data, count, MPI_DOUBLE, recv_data, count, MPI_INT, root, this->comm);
};
