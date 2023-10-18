#include <mutex>
#include <memory>
#include "mpi.h"

class FaradioMPI {
public:
    FaradioMPI();
    ~FaradioMPI();
    void init(int argc, char** argv);
    void finalize();

    static std::shared_ptr<FaradioMPI> getFaradioMPI();

    void send(double* data, int count, int dest, int tag);
    void recv(double* data, int count, int source, int tag);
    void reduce(double* send_data, double* recv_data, int count, int root);
    void reduce(void* send_data, void* recv_data, int count, int root);
    void gather(double* send_data, double* recv_data, int count, int root);

    int rank() const { return _rank; }
    int size() const { return _size; }
protected:
    MPI_Comm comm;
    int _rank;
    int _size;
};

static std::once_flag singletonFlag;
static std::shared_ptr<FaradioMPI> instance;