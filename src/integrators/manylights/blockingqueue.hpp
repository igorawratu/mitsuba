#ifndef BLOCKINGQUEUE_H_
#define BLOCKINGQUEUE_H_

#include <queue>
#include <mutex>
#include <condition_variable>
#include <tuple>

template <class T, class Container = std::deque<T>>
class BlockingQueue{
public:
    BlockingQueue(std::uint32_t max_size = 0) : max_size_(max_size), open_(true){}

    ~BlockingQueue(){}
    
    std::uint32_t size(){
        return queue_.size();
    }

    void close(){
        {
            std::unique_lock<std::mutex> lock(mutex_);
            open_ = false;
        }
        pull_cond_.notify_all();
    }

    void push(T&& item){
        {
            std::unique_lock<std::mutex> lock(mutex_);
            
            if(max_size_ > 0){
                push_cond_.wait(lock, [this](){return this->queue_.size() < max_size_;});
            }

            queue_.push(item);
        }

        pull_cond_.notify_one();
    }

    void push(const T& item){
        {
            std::unique_lock<std::mutex> lock(mutex_);
            
            if(max_size_ > 0){
                push_cond_.wait(lock, [this](){return this->queue_.size() < max_size_;});
            }

            queue_.push(item);
        }

        pull_cond_.notify_one();
    }

    bool pop(T& item){
        {
            std::unique_lock<std::mutex> lock(mutex_);

            pull_cond_.wait(lock, [this](){return !this->queue_.empty() || !this->open_;});
            if(queue_.empty()){
                return false;
            }
            
            item = std::move(queue_.front());
            queue_.pop();
        }

        push_cond_.notify_one();
        return true;
    }

private:
    std::queue<T, Container> queue_;
    std::mutex mutex_;
    std::condition_variable push_cond_;
    std::condition_variable pull_cond_;
    std::uint32_t max_size_;
    bool open_;
};

#endif