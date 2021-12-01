/*
    Copyright © 2021 Luiz Fernando F. G. Valle

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

/**
 * @file myvector.hpp
 * @author Luiz Fernando F. G. Valle
 * @brief Implements a vector template.
 *
 * @copyright Copyright © 2021 Luiz Fernando F. G. Valle
 */

#ifndef MYVECTOR_HPP
#define MYVECTOR_HPP

#include <cmath>   // for std::ceil
#include <cstring> // for std::memcpy
#include <memory>
#include <type_traits>

/**
 * @brief Contains the entire Vector module
 */
namespace myvector
{

/**
 * @brief Concepts used by the Vector implementation
 */
namespace concepts
{

/**
 * @brief Specifies that a type allows << operator with std::ostream as lhs
 *
 * @tparam T type
 */
template <typename T>
concept OstreamPrintable = requires(T obj, std::ostream ostream)
{
    ostream << obj;
};

} // namespace concepts

/**
 * @brief A vector template
 *
 * @tparam T Stored type
 */
template <typename T, typename Allocator = std::allocator<T>>
class Vector
{
    static_assert(std::is_same<T, typename Allocator::value_type>::value,
                  "Allocator must have T as its value_type");

    //                      //
    // Private type aliases //
    //                      //

    /**
     * @brief Alias for std::alocator_traits of the chosen allocator type
     */
    using alloctr_ = std::allocator_traits<Allocator>;

    //                                 //
    // Private static member variables //
    //                                 //

    /**
     * @brief The allocator used by the Vector
     */
    static Allocator alloc_;

public:
    //                     //
    // Public type aliases //
    //                     //

    /**
     * @brief Defining value_type as STL containers do
     */
    using value_type = T;

    /**
     * @brief Defining reference as STL containers do
     */
    using reference = T &;

    /**
     * @brief Defining const_reference as STL containers do
     */
    using const_reference = const T &;

    /**
     * @brief Defining allocator_type as STL containers do
     */
    using allocator_type = Allocator;

    /**
     * @brief Defining iterator as STL containers do
     */
    using iterator = T *;

    /**
     * @brief Defining const_iterator as STL containers do
     */
    using const_iterator = const T *;

    /**
     * @brief Defining pointer as STL containers do
     */
    using pointer = typename std::allocator_traits<Allocator>::pointer;

    /**
     * @brief Defining const_pointer as STL containers do
     */
    using const_pointer =
      typename std::allocator_traits<Allocator>::const_pointer;

    //                     //
    // RAII implementation //
    //                     //

    /**
     * @brief Constructs an empty Vector
     */
    Vector() = default;

    /**
     * @brief Constructs a Vector with the specified capacity
     *
     * @param capacity Number of elements
     */
    explicit Vector(size_t capacity)
        : elems_{alloctr_::allocate(alloc_, capacity)}, capacity_{capacity}
    {
    }

    /**
     * @brief Constructs a Vector with a list of elements
     *
     * @param list List of elements
     */
    Vector(std::initializer_list<T> list)
        : elems_{alloctr_::allocate(alloc_, list.size())},
          capacity_{list.size()},
          nelems_{list.size()}
    {
        static_assert(
          std::is_copy_constructible_v<T> || std::is_trivially_copyable_v<T>,
          "Vector initializer list constructor needs T to be copy"
          " constructible or trivially copyable");

        if constexpr (std::is_trivially_copyable_v<T>)
            std::memcpy(elems_, std::data(list), list.size() * sizeof(T));
        else
            for (auto dest{elems_}; const T &elem : list)
                alloctr_::construct(alloc_, dest++, elem);
    }

    /**
     * @brief Constructs a Vector by copy
     *
     * @param vec Vector to be copied
     */
    Vector(const Vector &vec) : nelems_{vec.nelems_}, capacity_{vec.capacity_}
    {
        elems_ = alloctr_::allocate(alloc_, capacity_);

        if constexpr (std::is_trivially_copyable_v<T>)
            std::memcpy(elems_, vec.elems_, nelems_ * sizeof(T));
        else
            for (auto dest{elems_}; auto &elem : vec)
                alloctr_::construct(alloc_, dest++, elem);
    }

    /**
     * @brief Constructs a Vector by move
     *
     * @param vec Vector to be moved from
     */
    Vector(Vector &&vec)
        : elems_{vec.elems_}, nelems_{vec.nelems_}, capacity_{vec.capacity_}
    {
        vec.elems_    = nullptr;
        vec.nelems_   = 0;
        vec.capacity_ = 0;
    }

    /**
     * @brief Destroys the Vector object
     */
    ~Vector() { DeallocateElems_(); }

    //                        //
    // Manipulation functions //
    //                        //

    /**
     * @brief Inserts an element into the end of the Vector by copy, increasing
     *          the Vector's size by 1.
     *
     * @param elem Element
     */
    void PushBack(const T &elem) requires std::copy_constructible<T>
    {
        GrowByFactorIfNeeded_();

        alloctr_::construct(alloc_, elems_ + nelems_++, elem);
    }

    /**
     * @brief Inserts an element into the end of the Vector by move, increasing
     *          the Vector's size by 1.
     *
     * @param elem Element
     */
    void PushBack(T &&elem) requires std::move_constructible<T>
    {
        GrowByFactorIfNeeded_();

        alloctr_::construct(alloc_, elems_ + nelems_++, std::move(elem));
    }

    /**
     * @brief Appends another Vector to the Vector.
     *
     * @param vector Vector to be copied from
     */
    void PushBack(const Vector<T> &vector) requires std::copy_constructible<T>
    {
        size_t total_nelems{nelems_ + vector.Size()};

        if (capacity_ < total_nelems)
            ResizeCapacity_(total_nelems);

        for (size_t i{nelems_}; auto &elem : vector)
            alloctr_::construct(alloc_, elems_ + i++, elem);

        nelems_ = total_nelems;
    }

    /**
     * @brief Appends another Vector to the Vector (by move).
     *
     * @param vector Vector to be moved from
     */
    void PushBack(Vector<T> &&vector) requires std::move_constructible<T>
    {
        size_t total_nelems{nelems_ + vector.Size()};

        if (capacity_ < total_nelems)
            ResizeCapacity_(total_nelems);

        for (size_t i{nelems_}; auto &elem : vector)
            alloctr_::construct(alloc_, elems_ + i++, std::move(elem));

        nelems_ = total_nelems;
    }

    /**
     * @brief Constructs and appends an element of type T to the Vector.
     *
     * @tparam Args Constructor arguments types
     * @param args Constructor arguments
     */
    template <typename... Args>
    void EmplaceBack(Args &&...args)
    {
        PushBack(T(std::forward<Args>(args)...));
    }

    /**
     * @brief Removes the last element from the Vector.
     *        This function should not be called with an empty Vector, checking
     *          is the caller's responsibility.
     */
    void PopBack()
    {
        --nelems_;

        if constexpr (!std::is_trivially_destructible<T>::value)
            alloctr_::destroy(alloc_, elems_ + nelems_);
    }

    /**
     * @brief Removes all elements from the Vector and calls their destructors
     *          if they are not trivially destructible.
     *        Does not free allocated memory.
     */
    void Clear()
    {
        if constexpr (!std::is_trivially_destructible<T>::value)
            DestroyElems_();

        nelems_ = 0;
    }

    //                                        //
    // Public capacity manipulation functions //
    //                                        //

    /**
     * @brief Reserves enough memory to store sz elements in the Vector, or
     *          does nothing if enough memory is already reserved.
     *        This function does not remove or modify any element already
     *          stored.
     *
     * @param sz Desired number of elements to be reserved in memory
     */
    void Reserve(size_t sz)
    {
        if (sz > capacity_)
            ResizeCapacity_(sz);
    }

    /**
     * @brief Shrinks the capacity to fit exactly the number of elements in the
     *          Vector. Does nothing if capacity_ is already at the minimum.
     */
    void ShrinkToFit()
    {
        if (capacity_ != nelems_)
            ResizeCapacity_(nelems_);
    }

    //                       //
    // Information functions //
    //                       //

    /**
     * @brief Returns the number of elements.
     *
     * @return Number of elements
     */
    size_t Size() const { return nelems_; }

    /**
     * @brief Checks if the Vector is empty
     *
     * @return true Vector is empty
     * @return false Vector is not empty
     */
    bool IsEmpty() const { return !nelems_; }

    /**
     * @brief Returns the Vector's current total capacity.
     *
     * @return Current capacity
     */
    size_t Capacity() const { return capacity_; }

    /**
     * @brief Returns the Vector's growth factor.
     *        Deprecated because this function will likely be removed.
     *
     * @return Vector's growth factor
     */
    [[deprecated]] float GrowthFactor() const { return GrowthFactor_; }

    //                  //
    // Setter functions //
    //                  //

    /**
     * @brief Change the Vector's growth factor.
     *        Does nothing if argument is <= 1.
     *        Deprecated because this function will likely be removed.
     *
     * @param factor New growth factor
     */
    [[deprecated]] void SetGrowthFactor(float factor)
    {
        if (factor > 1.f)
            GrowthFactor_ = factor;
    }

    //                    //
    // Iterator functions //
    //                    //

    /**
     * @brief Return an iterator to the first element.
     *
     * @return Pointer to the first element
     */
    pointer begin() const { return elems_; }

    /**
     * @brief Return an iterator to constant first element.
     *
     * @return Pointer to the constant first element
     */
    const_pointer cbegin() const { return elems_; }

    /**
     * @brief Return an iterator to one-past-last element.
     *
     * @return Pointer to the one-past-last element
     */
    pointer end() const { return elems_ + nelems_; }

    /**
     * @brief Return an iterator to constant one-past-last element.
     *
     * @return Pointer to the constant one-past-last element
     */
    const_pointer cend() const { return elems_ + nelems_; }

    //           //
    // Operators //
    //           //

    /*
        I want to reduce the following code duplication, but I'm not
          knowledgeable enough to do it yet.
        Attempts resulted in completely bizarre code.
    */

    /**
     * @brief Copies a Vector to another.
     *
     * @param rhs Vector to be copied from
     * @return Reference to left-hand Vector
     */
    Vector<T> &
    operator=(const Vector<T> &rhs) requires std::copy_constructible<T>
    {
        if (capacity_ >= rhs.Size())
        { // We can reuse the current Vector's memory
            size_t rhs_sz(rhs.Size());
            size_t min_sz{std::min(Size(), rhs_sz)};
            size_t elem_i{0};

            // For all elements already in the original Vector
            while (elem_i < min_sz)
            {
                if constexpr (std::is_copy_assignable<T>::value)
                    // Assign to them
                    elems_[elem_i] = rhs[elem_i];
                else // Write over them
                {
                    if constexpr (!std::is_trivially_destructible_v<T>)
                        // Need to destroy them first
                        alloctr_::destroy(alloc_, elems_ + elem_i);
                    alloctr_::construct(alloc_, elems_ + elem_i, rhs[elem_i]);
                }

                ++elem_i;
            }

            // Remaining elements are trash and don't need any treatment,
            //   so we just construct over them.
            while (elem_i < rhs_sz)
            {
                alloctr_::construct(alloc_, elems_ + elem_i, rhs[elem_i]);
                ++elem_i;
            }
        }
        else // Current memory is too small, needs reallocation
        {
            DeallocateElems_();
            ResizeCapacity_(rhs.Size());

            for (size_t i{nelems_}; const auto &elem : rhs)
                alloctr_::construct(alloc_, elems_ + i++, elem);
        }

        nelems_ = rhs.nelems_;

        return *this;
    }

    /**
     * @brief Moves a Vector to another.
     *
     * @param rhs Vector to be moved from
     * @return Reference to left-hand Vector
     */
    Vector<T> &operator=(Vector<T> &&rhs) requires std::move_constructible<T>
    {
        if (capacity_ >= rhs.Size())
        { // We can reuse the current Vector's memory
            size_t rhs_sz(rhs.Size());
            size_t min_sz{std::min(Size(), rhs_sz)};
            size_t elem_i{0};

            // For all elements already in the original Vector
            while (elem_i < min_sz)
            {
                if constexpr (std::is_move_assignable<T>::value)
                    // Assign to them
                    elems_[elem_i] = std::move(rhs[elem_i]);
                else // Write over them
                {
                    if constexpr (!std::is_trivially_destructible_v<T>)
                        // Need to destroy them first
                        alloctr_::destroy(alloc_, elems_ + elem_i);
                    alloctr_::construct(
                      alloc_, elems_ + elem_i, std::move(rhs[elem_i]));
                }

                ++elem_i;
            }

            // Remaining elements are trash and don't need any treatment,
            //   so we just construct over them.
            while (elem_i < rhs_sz)
            {
                alloctr_::construct(
                  alloc_, elems_ + elem_i, std::move(rhs[elem_i]));
                ++elem_i;
            }
        }
        else // Current memory is too small, needs reallocation
        {
            DeallocateElems_();
            ResizeCapacity_(rhs.Size());

            for (size_t i{nelems_}; auto &elem : rhs)
                alloctr_::construct(alloc_, elems_ + i++, std::move(elem));
        }

        nelems_ = rhs.nelems_;

        rhs.Clear();

        return *this;
    }

    /**
     * @brief Returns a reference to a specified element
     *
     * @param i Index of the element
     * @return Reference to the element
     */
    reference operator[](size_t i) { return elems_[i]; }

    /**
     * @brief Returns a const reference to a specified element
     *
     * @param i Index of the element
     * @return Const reference to the element
     */
    const_reference &operator[](size_t i) const { return elems_[i]; }

    /**
     * @brief Spaceship operator on two Vector's number of elements
     *
     * @param rhs Vector to compare with
     * @return Spaceship operation result
     */
    auto operator<=>(const Vector<T> &rhs) const
    {
        return nelems_ <=> rhs.nelems_;
    }

    /*
        I'm still thinking whether or not using templates for the operators
          is better than hard coding them, but it results in less lines of code
          so why not.
    */

    /**
     * @brief Appends to the Vector.
     *        Wraps around PushBack().
     *
     * @tparam RhsT Type of the argument
     * @param elem Argument given to PushBack
     * @return Reference to the Vector
     */
    template <typename RhsT>
    Vector<T> &operator+=(const RhsT &elem)
    {
        PushBack(elem);

        return *this;
    }

    /**
     * @brief Appends the argument to a copy of the vector and returns the copy.
     *        Wraps around operator+=().
     *
     * @tparam RhsT Type of the argument
     * @param elem Argument given to operator+=()
     * @return Vector<T> Resulting Vector
     */
    template <typename RhsT>
    Vector<T> operator+(const RhsT &elem) const
    {
        Vector<T> ret{*this};

        ret += elem;

        return ret;
    }

private:
    //                                   //
    // Private element removal functions //
    //                                   //

    /**
     * @brief Destroys all elements in the Vector.
     *        Disabled if T is trivially destructible.
     *
     * @tparam Ret Return type, shouldn't be manually specified.
     */
    template <typename Ret = void>
    std::enable_if_t<!std::is_trivially_destructible<T>::value,
                     Ret> constexpr DestroyElems_()
    {
        for (auto *end_ptr{elems_ + nelems_}, *elem_ptr{elems_};
             elem_ptr < end_ptr;
             ++elem_ptr)
            alloctr_::destroy(alloc_, elem_ptr);

        nelems_ = 0;
    }

    /**
     * @brief Deallocates the Vector's representation.
     *        Also calls destructor for each element if not trivially
     *        destructible.
     */
    void DeallocateElems_()
    {
        if constexpr (!std::is_trivially_destructible<T>::value)
            DestroyElems_();
        else
            nelems_ = 0;

        alloctr_::deallocate(alloc_, elems_, capacity_);
        capacity_ = 0;
    }

    //                                 //
    // Capacity manipulation functions //
    //                                 //

    /**
     * @brief Resizes the capacity to the desired size.
     *        The Vector will keep all elements that fit in the new size,
     *          maintaining their order.
     *
     * @param sz Desired size
     */
    void ResizeCapacity_(size_t sz) requires std::move_constructible<T>
    {
        if constexpr (std::copy_constructible<T>)
        {
            Vector<T> old{*this};

            capacity_ = sz;
            elems_    = alloctr_::allocate(alloc_, capacity_);
            nelems_   = std::min(old.nelems_, capacity_);

            for (size_t i{0}; auto &elem : old)
                alloctr_::construct(alloc_, elems_ + i++, elem);
        }
        else
        {
            Vector<T> old{std::move(*this)};

            capacity_ = sz;
            elems_    = alloctr_::allocate(alloc_, capacity_);
            nelems_   = std::min(old.nelems_, capacity_);

            for (size_t i{0}; auto &elem : old)
                alloctr_::construct(alloc_, elems_ + i++, std::move(elem));
        }
    }

    /**
     * @brief Grows the Vector's capacity multiplying it by GrowthFactor_.
     *        Wraps around ResizeCapacity_.
     */
    void GrowByFactor_()
    {
        ResizeCapacity_(capacity_ ? std::ceil(capacity_ * double{GrowthFactor_})
                                  : 1);
    }

    /**
     * @brief Grows the Vector's capacity if full.
     *        Wraps around GrowByFactor_.
     */
    void GrowByFactorIfNeeded_()
    {
        if (IsFullCapacity_())
            GrowByFactor_();
    }

    //                               //
    // Private information functions //
    //                               //

    /**
     * @brief Checks if capacity is full.
     *
     * @return true Vector is using all its capacity
     * @return false Vector is not using all its capacity
     */
    bool IsFullCapacity_() const { return nelems_ >= capacity_; }

    //                          //
    // Private member variables //
    //                          //

    /**
     * @brief Pointer to the array of elements
     */
    pointer elems_{nullptr};

    /**
     * @brief Number of elements in the vector
     */
    size_t nelems_{0};

    /**
     * @brief Number of elements allocated or preallocated in the vector
     */
    size_t capacity_{0};

    /*
        Unsure if GrowthFactor_ should be static. I think having it static
          is probably best, but I'll see from experience.
        Its setter and getters are marked as deprecated becaue of their possible
          future removal.
    */

    /**
     * @brief The allocated size will be multiplied by
     *          this value when necessary.
     *        Must be greater than 1
     */
    float GrowthFactor_{1.5}; // Or, more precisely, 1.61803 (phi)
};

//                               //
// Static members initialization //
//                               //

template <typename T, typename Allocator>
Allocator Vector<T, Allocator>::alloc_{};

//                          //
// Other function overloads //
//                          //

/**
 * @brief Outputs a Vector to an std::ostream.
 *
 * @tparam T Vector's value type
 * @tparam Allocator Vector's allocator type
 * @param lhs Ostream target
 * @param rhs Vector source
 * @return std::ostream& Reference to lhs std::ostream
 */
template <typename T, typename Allocator>
std::ostream &operator<<(std::ostream &lhs, Vector<T, Allocator> &rhs)
{
    static_assert(concepts::OstreamPrintable<T>,
                  "Can't output Vector to ostream with operator << "
                  "because the underlying type has no such operation.");

    // Reduce useless errors by ignoring the following code
    //   if the above assertion fails
    if constexpr (concepts::OstreamPrintable<T>)
    {
        lhs << '[';

        if (rhs.Size())
            lhs << rhs[0];

        for (size_t elem_i{1}; elem_i < rhs.Size(); ++elem_i)
            lhs << ", " << rhs[elem_i];

        lhs << ']';
    }

    return lhs;
}

} // namespace myvector

#endif // #ifndef MYVECTOR_HPP