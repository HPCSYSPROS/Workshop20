/*
All modification made by Intel Corporation: © 2016 Intel Corporation

All contributions by the University of California:
Copyright (c) 2014, 2015, The Regents of the University of California (Regents)
All rights reserved.

All other contributions:
Copyright (c) 2014, 2015, the respective contributors
All rights reserved.
For the list of contributors go to https://github.com/BVLC/caffe/blob/master/CONTRIBUTORS.md


Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <string>
#include <vector>

#include "caffe/sgd_solvers.hpp"
#include "caffe/util/hdf5.hpp"
#include "caffe/util/io.hpp"
#include "caffe/util/upgrade_proto.hpp"

namespace caffe {
template <typename Dtype>
Dtype SGDSolver<Dtype>::GetWarmUpLR(int cur_iter, int warmup_iter, Dtype warmup_start_lr) {
  if (cur_iter < 0) {
    cur_iter = 0;
  }
  return (cur_iter * this->param_.base_lr() +
          (warmup_iter - cur_iter) * warmup_start_lr) / warmup_iter;
}
// Return the current learning rate. The currently implemented learning rate
// policies are as follows:
//    - fixed: always return base_lr.
//    - step: return base_lr * gamma ^ (floor(iter / step))
//    - exp: return base_lr * gamma ^ iter
//    - inv: return base_lr * (1 + gamma * iter) ^ (- power)
//    - multistep: similar to step but it allows non uniform steps defined by
//      stepvalue
//    - poly: the effective learning rate follows a polynomial decay, to be
//      zero by the max_iter. return base_lr (1 - iter/max_iter) ^ (power)
//    - sigmoid: the effective learning rate follows a sigmod decay
//      return base_lr ( 1/(1 + exp(-gamma * (iter - stepsize))))
//
// where base_lr, max_iter, gamma, step, stepvalue and power are defined
// in the solver parameter protocol buffer, and iter is the current iteration.
template <typename Dtype>
Dtype SGDSolver<Dtype>::GetLearningRate() {
  Dtype rate;
  const string& lr_policy = this->param_.lr_policy();


  if (this->param_.warmup_iter() > 0 &&
      this->iter_ < this->param_.warmup_iter()) {
    rate = GetWarmUpLR(this->iter_, this->param_.warmup_iter(),
                       this->param_.warmup_start_lr());
  } else if (lr_policy == "fixed") {
    rate = this->param_.base_lr();
  } else if (lr_policy == "step") {
    this->current_step_ = this->iter_ / this->param_.stepsize();
    rate = this->param_.base_lr() *
        pow(this->param_.gamma(), this->current_step_);
  } else if (lr_policy == "exp") {
    rate = this->param_.base_lr() * pow(this->param_.gamma(), this->iter_);
  } else if (lr_policy == "inv") {
    rate = this->param_.base_lr() *
        pow(Dtype(1) + this->param_.gamma() * this->iter_,
            - this->param_.power());
  } else if (lr_policy == "multistep") {
    if (this->current_step_ < this->param_.stepvalue_size() &&
          this->iter_ >= this->param_.stepvalue(this->current_step_)) {
      this->current_step_++;
      LOG(INFO) << "MultiStep Status: Iteration " <<
      this->iter_ << ", step = " << this->current_step_;
    }
    rate = this->param_.base_lr() *
        pow(this->param_.gamma(), this->current_step_);
  } else if (lr_policy == "poly") {
    rate = this->param_.base_lr() * pow(Dtype(1.) -
        (Dtype(this->iter_) / Dtype(this->param_.max_iter())),
        this->param_.power());
  } else if (lr_policy == "sigmoid") {
    rate = this->param_.base_lr() * (Dtype(1.) /
        (Dtype(1.) + exp(-this->param_.gamma() * (Dtype(this->iter_) -
          Dtype(this->param_.stepsize())))));
  } else if (lr_policy == "plateau") {
    // Update minimum loss if needed
    if (this->smoothed_loss_ < this->minimum_loss_) {
      this->minimum_loss_ = this->smoothed_loss_;
      this->iter_last_event_ = this->iter_;
    }

    // If sufficient iters have passed after the last event, then lower LR
    // An event is defined an update of minimum loss or LR
    if (this->current_step_ < this->param_.plateau_winsize_size()) {
      int iter_next_update = this->iter_last_event_
            + this->param_.plateau_winsize(this->current_step_);

      if (this->iter_ >= iter_next_update) {
        this->current_step_++;
        this->iter_last_event_ = this->iter_;
        LOG(INFO) << "Plateau Status: Iteration " << this->iter_
                  << ", step = " << this->current_step_;
      }
    }

    if (this->param_.display() && this->iter_ % this->param_.display() == 0
        && this->iter_last_event_ > (this->iter_ - this->param_.display())) {
      LOG(INFO) << "Plateau Status: Iteration " << this->iter_
                << ", current minimum_loss = " << this->minimum_loss_;
    }

    rate = this->param_.base_lr() *
        pow(this->param_.gamma(), this->current_step_);
  } else {
    LOG(FATAL) << "Unknown learning rate policy: " << lr_policy;
  }
  return rate;
}

template <typename Dtype>
void SGDSolver<Dtype>::PreSolve() {
  // Initialize the history
  const vector<Blob<Dtype>*>& net_params = this->net_->learnable_params();
  history_.clear();
  update_.clear();
  temp_.clear();
  for (int i = 0; i < net_params.size(); ++i) {
    const vector<int>& shape = net_params[i]->shape();

    // TODO: allocate these buffers taking into account owned_count to reduce memory footprint
    history_.push_back(shared_ptr<Blob<Dtype> >(new Blob<Dtype>(shape)));
    update_.push_back(shared_ptr<Blob<Dtype> >(new Blob<Dtype>(shape)));
    temp_.push_back(shared_ptr<Blob<Dtype> >(new Blob<Dtype>(shape)));
  }

  this->minimum_loss_ = std::numeric_limits<float>::max();
}

template <typename Dtype>
void SGDSolver<Dtype>::ClipGradients() {
  const Dtype clip_gradients = this->param_.clip_gradients();
  if (clip_gradients < 0) { return; }
  const vector<Blob<Dtype>*>& net_params = this->net_->learnable_params();
  Dtype sumsq_diff = 0;
  for (int i = 0; i < net_params.size(); ++i) {
    sumsq_diff += net_params[i]->sumsq_diff();
  }
  const Dtype l2norm_diff = std::sqrt(sumsq_diff);
  if (l2norm_diff > clip_gradients) {
    Dtype scale_factor = clip_gradients / l2norm_diff;
    LOG(INFO) << "Gradient clipping: scaling down gradients (L2 norm "
        << l2norm_diff << " > " << clip_gradients << ") "
        << "by scale factor " << scale_factor;
    for (int i = 0; i < net_params.size(); ++i) {
      net_params[i]->scale_diff(scale_factor);
    }
  }
}

template <typename Dtype>
void SGDSolver<Dtype>::ApplyUpdate() {
  CHECK(Caffe::root_solver());
  Dtype rate = GetLearningRate();
  if (this->param_.display() && this->iter_ % this->param_.display() == 0) {
    LOG(INFO) << "Iteration " << this->iter_ << ", lr = " << rate;
  }
  ClipGradients();
  for (int param_id = 0; param_id < this->net_->learnable_params().size();
       ++param_id) {
    ApplyUpdate(param_id);
  }
}

template <typename Dtype>
void SGDSolver<Dtype>::ApplyUpdate(int param_id) {
  CHECK(Caffe::root_solver());
  Dtype rate = GetLearningRate();

  LOG_PARAM_BLOB(this->net_->learnable_params()[param_id], diff, param_id, "ApplyUpdate: raw delwt:");

  // If Learning rate for this learnable params is zero then skip
  // updating params
  if (this->net_->params_lr()[param_id] == 0) {
    return;
  }

  Normalize(param_id);
  LOG_PARAM_BLOB(this->net_->learnable_params()[param_id], diff, param_id, "ApplyUpdate: delwt after Normalize:");

  Regularize(param_id);
  LOG_PARAM_BLOB(this->net_->learnable_params()[param_id], diff, param_id, "ApplyUpdate: delwt after Regularize:");

  ComputeUpdateValue(param_id, rate);
  LOG_PARAM_BLOB(this->net_->learnable_params()[param_id], diff, param_id, "ApplyUpdate: wtinc:");

  LOG_PARAM_BLOB(this->net_->learnable_params()[param_id], data, param_id, "ApplyUpdate: weight before update:");

  this->net_->learnable_params()[param_id]->Update();

  LOG_PARAM_BLOB(this->net_->learnable_params()[param_id], data, param_id, "ApplyUpdate: weight after update:");
}

template <typename Dtype>
void SGDSolver<Dtype>::Normalize(int param_id) {

#ifdef USE_MLSL
  if ((this->param_.iter_size() == 1) && !mn::is_multinode()) {
    return;
  }
#else /* !USE_MLSL */
  if (this->param_.iter_size() == 1) { return; }
#endif /* USE_MLSL */

  // Scale gradient to counterbalance accumulation.
  const vector<Blob<Dtype>*>& net_params = this->net_->learnable_params();

#ifdef USE_MLSL
  const Dtype accum_normalization = Dtype(1.) / (this->param_.iter_size() * mn::get_nodes_count());
#else /* !USE_MLSL */
  const Dtype accum_normalization = Dtype(1.) / this->param_.iter_size();
#endif /* USE_MLSL */

  switch (Caffe::mode()) {
  case Caffe::CPU: {

    if (net_params[param_id]->prv_diff()
        && (net_params[param_id]->prv_diff_count()
            == net_params[param_id]->count())) {

        caffe_scal(net_params[param_id]->count(), accum_normalization,
            net_params[param_id]->mutable_prv_diff());
    }
    else {
        caffe_scal(net_params[param_id]->count(), accum_normalization,
            net_params[param_id]->mutable_cpu_diff());
    }

    break;
  }
  case Caffe::GPU: {
#ifndef CPU_ONLY
    caffe_gpu_scal(net_params[param_id]->count(), accum_normalization,
        net_params[param_id]->mutable_gpu_diff());
#else
    NO_GPU;
#endif
    break;
  }
  default:
    LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
  }
}

template <typename Dtype>
void SGDSolver<Dtype>::Regularize(int param_id) {
  const vector<Blob<Dtype>*>& net_params = this->net_->learnable_params();
  const vector<float>& net_params_weight_decay =
      this->net_->params_weight_decay();
  Dtype weight_decay = this->param_.weight_decay();
  string regularization_type = this->param_.regularization_type();
  Dtype local_decay = weight_decay * net_params_weight_decay[param_id];
  switch (Caffe::mode()) {
  case Caffe::CPU: {
    if (local_decay) {
      if (regularization_type == "L2") {
        // add weight decay
        if (net_params[param_id]->prv_data()
             && (net_params[param_id]->prv_data_count()
                 == net_params[param_id]->count())) {
          CHECK_EQ(true,
            net_params[param_id]->get_prv_data_descriptor()->layout_compare(
            net_params[param_id]->get_prv_diff_descriptor()));

          caffe_axpy(net_params[param_id]->count(),
                     local_decay,
                     net_params[param_id]->prv_data(),
                     net_params[param_id]->mutable_prv_diff());
        } else {
          caffe_axpy(net_params[param_id]->count(),
              local_decay,
              net_params[param_id]->cpu_data(),
              net_params[param_id]->mutable_cpu_diff());
        }
      } else if (regularization_type == "L1") {
        caffe_cpu_sign(net_params[param_id]->count(),
            net_params[param_id]->cpu_data(),
            temp_[param_id]->mutable_cpu_data());
        caffe_axpy(net_params[param_id]->count(),
            local_decay,
            temp_[param_id]->cpu_data(),
            net_params[param_id]->mutable_cpu_diff());
      } else {
        LOG(FATAL) << "Unknown regularization type: " << regularization_type;
      }
    }
    break;
  }
  case Caffe::GPU: {
#ifndef CPU_ONLY
    if (local_decay) {
      if (regularization_type == "L2") {
        // add weight decay
        caffe_gpu_axpy(net_params[param_id]->count(),
            local_decay,
            net_params[param_id]->gpu_data(),
            net_params[param_id]->mutable_gpu_diff());
      } else if (regularization_type == "L1") {
        caffe_gpu_sign(net_params[param_id]->count(),
            net_params[param_id]->gpu_data(),
            temp_[param_id]->mutable_gpu_data());
        caffe_gpu_axpy(net_params[param_id]->count(),
            local_decay,
            temp_[param_id]->gpu_data(),
            net_params[param_id]->mutable_gpu_diff());
      } else {
        LOG(FATAL) << "Unknown regularization type: " << regularization_type;
      }
    }
#else
    NO_GPU;
#endif
    break;
  }
  default:
    LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
  }
}

#ifndef CPU_ONLY
template <typename Dtype>
void sgd_update_gpu(int N, Dtype* g, Dtype* h, Dtype momentum,
    Dtype local_rate);
#endif

template <typename Dtype>
void SGDSolver<Dtype>::ComputeUpdateValue(int param_id, Dtype rate) {
  const vector<Blob<Dtype>*>& net_params = this->net_->learnable_params();
  const vector<float>& net_params_lr = this->net_->params_lr();
  Dtype momentum = this->param_.momentum();
  Dtype local_rate = rate * net_params_lr[param_id];

  if (this->param_.warmup_iter() > 0 &&
      this->iter_ < this->param_.warmup_iter()) {
    // Momentum correction during warmup stage
    Dtype prev_rate = GetWarmUpLR(this->iter_ - 1, this->param_.warmup_iter(),
                                  this->param_.warmup_start_lr());
    momentum = momentum * (rate / prev_rate);
  }
  // Compute the update to history, then copy it to the parameter diff.
  switch (Caffe::mode()) {
  case Caffe::CPU: {
    if (net_params[param_id]->prv_diff()
        && (net_params[param_id]->prv_diff_count()
            == net_params[param_id]->count())) {

      caffe_cpu_axpby(net_params[param_id]->count(), local_rate,
                      net_params[param_id]->prv_diff(), momentum,
                      history_[param_id]->mutable_cpu_data());

      caffe_copy(net_params[param_id]->count(),
                 history_[param_id]->cpu_data(),
                 net_params[param_id]->mutable_prv_diff());
    } else {
      caffe_cpu_axpby(net_params[param_id]->count(), local_rate,
                     net_params[param_id]->cpu_diff(), momentum,
                     history_[param_id]->mutable_cpu_data());

      caffe_copy(net_params[param_id]->count(),
                 history_[param_id]->cpu_data(),
                 net_params[param_id]->mutable_cpu_diff());
    }
    break;
  }
  case Caffe::GPU: {
#ifndef CPU_ONLY
    sgd_update_gpu(net_params[param_id]->count(),
        net_params[param_id]->mutable_gpu_diff(),
        history_[param_id]->mutable_gpu_data(),
        momentum, local_rate);
#else
    NO_GPU;
#endif
    break;
  }
  default:
    LOG(FATAL) << "Unknown caffe mode: " << Caffe::mode();
  }
}

template <typename Dtype>
void SGDSolver<Dtype>::SnapshotSolverState(const string& model_filename) {
  switch (this->param_.snapshot_format()) {
    case caffe::SolverParameter_SnapshotFormat_BINARYPROTO:
      SnapshotSolverStateToBinaryProto(model_filename);
      break;
    case caffe::SolverParameter_SnapshotFormat_HDF5:
      SnapshotSolverStateToHDF5(model_filename);
      break;
    default:
      LOG(FATAL) << "Unsupported snapshot format.";
  }
}

template <typename Dtype>
void SGDSolver<Dtype>::SnapshotSolverStateToBinaryProto(
    const string& model_filename) {
  SolverState state;
  state.set_iter(this->iter_);
  state.set_learned_net(model_filename);
  state.set_current_step(this->current_step_);
  state.set_iter_last_event(this->iter_last_event_);
  state.set_minimum_loss(this->minimum_loss_);
  state.clear_history();
  for (int i = 0; i < history_.size(); ++i) {
    // Add history
    BlobProto* history_blob = state.add_history();
    history_[i]->ToProto(history_blob);
  }
  string snapshot_filename = Solver<Dtype>::SnapshotFilename(".solverstate");
#ifdef USE_MLSL
  if (mn::is_root()) {
#endif
  LOG(INFO)
    << "Snapshotting solver state to binary proto file " << snapshot_filename;
  WriteProtoToBinaryFile(state, snapshot_filename.c_str());
#ifdef USE_MLSL
  }
#endif
}

template <typename Dtype>
void SGDSolver<Dtype>::SnapshotSolverStateToHDF5(
    const string& model_filename) {
  string snapshot_filename =
      Solver<Dtype>::SnapshotFilename(".solverstate.h5");
  LOG(INFO) << "Snapshotting solver state to HDF5 file " << snapshot_filename;
  hid_t file_hid = H5Fcreate(snapshot_filename.c_str(), H5F_ACC_TRUNC,
      H5P_DEFAULT, H5P_DEFAULT);
  CHECK_GE(file_hid, 0)
      << "Couldn't open " << snapshot_filename << " to save solver state.";
  hdf5_save_int(file_hid, "iter", this->iter_);
  hdf5_save_string(file_hid, "learned_net", model_filename);
  hdf5_save_int(file_hid, "current_step", this->current_step_);
  hdf5_save_int(file_hid, "iter_last_event", this->iter_last_event_);
  hdf5_save_float<Dtype>(file_hid, "minimum_loss", this->minimum_loss_);
  hid_t history_hid = H5Gcreate2(file_hid, "history", H5P_DEFAULT, H5P_DEFAULT,
      H5P_DEFAULT);
  CHECK_GE(history_hid, 0)
      << "Error saving solver state to " << snapshot_filename << ".";
  for (int i = 0; i < history_.size(); ++i) {
    ostringstream oss;
    oss << i;
    hdf5_save_nd_dataset<Dtype>(history_hid, oss.str(), *history_[i]);
  }
  H5Gclose(history_hid);
  H5Fclose(file_hid);
}

template <typename Dtype>
void SGDSolver<Dtype>::RestoreSolverStateFromBinaryProto(
    const string& state_file) {
  SolverState state;
  ReadProtoFromBinaryFile(state_file, &state);
  this->iter_ = state.iter();
  if (state.has_learned_net()) {
    NetParameter net_param;
    ReadNetParamsFromBinaryFileOrDie(state.learned_net().c_str(), &net_param);
    this->net_->CopyTrainedLayersFrom(net_param);
  }
  this->current_step_ = state.current_step();
  this->iter_last_event_ = state.iter_last_event();
  this->minimum_loss_ = state.minimum_loss();
  CHECK_EQ(state.history_size(), history_.size())
      << "Incorrect length of history blobs.";
  LOG(INFO) << "SGDSolver: restoring history";
  for (int i = 0; i < history_.size(); ++i) {
    history_[i]->FromProto(state.history(i));
  }
}

template <typename Dtype>
void SGDSolver<Dtype>::RestoreSolverStateFromHDF5(const string& state_file) {
  hid_t file_hid = H5Fopen(state_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  CHECK_GE(file_hid, 0) << "Couldn't open solver state file " << state_file;
  this->iter_ = hdf5_load_int(file_hid, "iter");
  if (H5LTfind_dataset(file_hid, "learned_net")) {
    string learned_net = hdf5_load_string(file_hid, "learned_net");
    this->net_->CopyTrainedLayersFrom(learned_net);
  }
  this->current_step_ = hdf5_load_int(file_hid, "current_step");
  this->iter_last_event_ = hdf5_load_int(file_hid, "iter_last_event");
  this->minimum_loss_ = hdf5_load_float<Dtype>(file_hid, "minimum_loss");
  hid_t history_hid = H5Gopen2(file_hid, "history", H5P_DEFAULT);
  CHECK_GE(history_hid, 0) << "Error reading history from " << state_file;
  int state_history_size = hdf5_get_num_links(history_hid);
  CHECK_EQ(state_history_size, history_.size())
      << "Incorrect length of history blobs.";
  for (int i = 0; i < history_.size(); ++i) {
    ostringstream oss;
    oss << i;
    hdf5_load_nd_dataset<Dtype>(history_hid, oss.str().c_str(), 0,
                                kMaxBlobAxes, history_[i].get());
  }
  H5Gclose(history_hid);
  H5Fclose(file_hid);
}

INSTANTIATE_CLASS(SGDSolver);
REGISTER_SOLVER_CLASS(SGD);

}  // namespace caffe
