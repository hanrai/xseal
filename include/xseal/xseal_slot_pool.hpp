// SPDX-License-Identifier: Apache-2.0
// Copyright (c) 2026 hanrai. All Rights Reserved.

#pragma once

#include <atomic>
#include <cstdint>
#include <vector>

namespace xseal {

template <uint32_t NumSlots = 768> class XsealSlotPool;

/**
 * XsealSlot: RAII handle for a buffer slot.
 * When the last copy of this handle is destroyed, the slot is automatically
 * returned to the pool's free queue.
 */
template <uint32_t NumSlots = 768> class XsealSlot {
public:
  XsealSlot() : slot_id(-1), pool(nullptr) {}
  XsealSlot(int32_t id, XsealSlotPool<NumSlots> *p);

  ~XsealSlot();

  // Copying increments the reference count
  XsealSlot(const XsealSlot &other);
  XsealSlot &operator=(const XsealSlot &other);

  // Moving transfers ownership without touching refcount
  XsealSlot(XsealSlot &&other) noexcept;
  XsealSlot &operator=(XsealSlot &&other) noexcept;

  int32_t id() const { return slot_id; }
  bool valid() const { return slot_id != -1 && pool != nullptr; }

  void reset();

private:
  int32_t slot_id;
  XsealSlotPool<NumSlots> *pool;
};

/**
 * XsealSlotPool: Manages buffer slots with cache-line aligned counters.
 */
template <uint32_t NumSlots> class XsealSlotPool {
public:
  static constexpr uint32_t NUM_SLOTS = NumSlots;

  XsealSlotPool() {
    for (uint32_t i = 0; i < NUM_SLOTS; ++i) {
      counters[i].val.store(0, std::memory_order_relaxed);
      free_queue[i].store(i, std::memory_order_relaxed);
    }
    q_head.store(0, std::memory_order_relaxed);
    q_tail.store(NUM_SLOTS, std::memory_order_relaxed);
  }

  /**
   * Acquire a slot from the free queue.
   * Returns a valid XsealSlot on success, or an invalid one if pool is empty.
   */
  XsealSlot<NumSlots> acquire() {
    uint32_t h = q_head.load(std::memory_order_acquire);
    while (h != q_tail.load(std::memory_order_acquire)) {
      if (q_head.compare_exchange_weak(h, h + 1, std::memory_order_release)) {
        int32_t id = free_queue[h % NUM_SLOTS].load(std::memory_order_relaxed);
        // Initial lease for the reader/parser
        counters[id].val.store(1, std::memory_order_relaxed);
        return XsealSlot<NumSlots>(id, this);
      }
    }
    return XsealSlot<NumSlots>(); // Empty
  }

  // Internal API for XsealSlot
  void add_ref(int32_t id) {
    if (id >= 0 && id < (int32_t)NUM_SLOTS) {
      counters[id].val.fetch_add(1, std::memory_order_relaxed);
    }
  }

  void release(int32_t id) {
    if (id >= 0 && id < (int32_t)NUM_SLOTS) {
      if (counters[id].val.fetch_sub(1, std::memory_order_acq_rel) == 1) {
        // Last one out, return to queue
        uint32_t t = q_tail.fetch_add(1, std::memory_order_acq_rel);
        free_queue[t % NUM_SLOTS].store(id, std::memory_order_release);
      }
    }
  }

private:
  struct alignas(64) AlignedCounter {
    std::atomic<int32_t> val;
  };

  AlignedCounter counters[NUM_SLOTS];

  // Simple lock-free MPSC-ish queue for indices
  std::atomic<int32_t> free_queue[NUM_SLOTS];
  std::atomic<uint32_t> q_head;
  std::atomic<uint32_t> q_tail;
};

// Inline implementations for XsealSlot
template <uint32_t NumSlots>
inline XsealSlot<NumSlots>::XsealSlot(int32_t id, XsealSlotPool<NumSlots> *p)
    : slot_id(id), pool(p) {}

template <uint32_t NumSlots> inline XsealSlot<NumSlots>::~XsealSlot() { reset(); }

template <uint32_t NumSlots>
inline XsealSlot<NumSlots>::XsealSlot(const XsealSlot<NumSlots> &other)
    : slot_id(other.slot_id), pool(other.pool) {
  if (valid())
    pool->add_ref(slot_id);
}

template <uint32_t NumSlots>
inline XsealSlot<NumSlots> &XsealSlot<NumSlots>::operator=(const XsealSlot<NumSlots> &other) {
  if (this != &other) {
    reset();
    slot_id = other.slot_id;
    pool = other.pool;
    if (valid())
      pool->add_ref(slot_id);
  }
  return *this;
}

template <uint32_t NumSlots>
inline XsealSlot<NumSlots>::XsealSlot(XsealSlot<NumSlots> &&other) noexcept
    : slot_id(other.slot_id), pool(other.pool) {
  other.slot_id = -1;
  other.pool = nullptr;
}

template <uint32_t NumSlots>
inline XsealSlot<NumSlots> &XsealSlot<NumSlots>::operator=(XsealSlot<NumSlots> &&other) noexcept {
  if (this != &other) {
    reset();
    slot_id = other.slot_id;
    pool = other.pool;
    other.slot_id = -1;
    other.pool = nullptr;
  }
  return *this;
}

template <uint32_t NumSlots> inline void XsealSlot<NumSlots>::reset() {
  if (valid()) {
    pool->release(slot_id);
    slot_id = -1;
    pool = nullptr;
  }
}

} // namespace xseal
