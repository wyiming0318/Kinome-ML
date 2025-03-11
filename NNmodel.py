import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import math, copy, time
from torch.autograd import Variable
import ndf


class SinusoidalPosEmb(torch.nn.Module):
	def __init__(self, dim):
		super().__init__()
		self.dim = dim

	def forward(self, x):
		device = x.device
		half_dim = self.dim // 2
		emb = math.log(10000) / (half_dim - 1)
		emb = torch.exp(torch.arange(half_dim, device=device) * -emb)
		emb = x[:, None] * emb[None, :]
		emb = torch.cat((emb.sin(), emb.cos()), dim=-1)
		return emb


class MemNN(torch.nn.Module):
	def __init__(self, dim_main, dim_joint, h_list_main, h_list_joint1, h_list_joint2, dropout, n_tree, tree_depth, tree_feature_rate):
		super().__init__()
		self.dim_main = dim_main
		self.dim_joint = dim_joint
		self.num_layers_main = len(h_list_main)
		self.num_layers_joint1 = len(h_list_joint1)
		self.num_layers_joint2 = len(h_list_joint2)
		self.dropout = dropout

		self.forest = ndf.Forest(n_tree=n_tree, tree_depth=tree_depth, n_in_feature=h_list_main[0],
								tree_feature_rate=tree_feature_rate, n_class=2,
								jointly_training=True)


		#MLP_aux1 and MLP_aux2
		self.layer_num_emb = SinusoidalPosEmb(h_list_joint1[-1])
		self.linear_layer_num_emb_scale = nn.Linear(h_list_joint1[-1], h_list_joint1[-1])
		self.linear_layer_num_emb_shift = nn.Linear(h_list_joint1[-1], h_list_joint1[-1])
		self.dropout_layer_num_emb_scale = nn.Dropout(p=dropout)
		self.dropout_layer_num_emb_shift = nn.Dropout(p=dropout)


		#MLP_main
		self.norm_main = nn.ModuleList([])
		self.linears_main = nn.ModuleList([])
		for i in range(self.num_layers_main):
			if i == 0:
				self.linears_main.append(nn.Linear(self.dim_main, h_list_main[i]))
			else:
				self.linears_main.append(nn.Linear(h_list_main[i-1], h_list_main[i]))
			self.norm_main.append(nn.LayerNorm(h_list_main[i]))

		self.dropout_main = nn.ModuleList([])
		for i in range(self.num_layers_main-1):
			self.dropout_main.append(nn.Dropout(p=dropout))


		#MLP_p1
		self.norm_joint1 = nn.ModuleList([])
		self.linears_joint1 = nn.ModuleList([])
		for i in range(self.num_layers_joint1):
			if i == 0:
				self.linears_joint1.append(nn.Linear(self.dim_joint, h_list_joint1[i]))
			else:
				self.linears_joint1.append(nn.Linear(h_list_joint1[i-1], h_list_joint1[i]))
			self.norm_joint1.append(nn.LayerNorm(h_list_joint1[i]))

		self.dropout_joint1 = nn.ModuleList([])
		for i in range(self.num_layers_joint1):
			self.dropout_joint1.append(nn.Dropout(p=dropout))


		#Shared part of MLP_p2 and MLP_p3
		self.norm_joint2 = nn.ModuleList([])
		self.linears_joint2 = nn.ModuleList([])
		for i in range(self.num_layers_joint2):
			if i == 0:
				self.linears_joint2.append(nn.Linear(h_list_joint1[-1], h_list_joint2[i]))
			else:
				self.linears_joint2.append(nn.Linear(h_list_joint2[i-1], h_list_joint2[i]))
			self.norm_joint2.append(nn.LayerNorm(h_list_joint2[i]))

		self.dropout_joint2 = nn.ModuleList([])
		for i in range(self.num_layers_joint2):
			self.dropout_joint2.append(nn.Dropout(p=dropout))


		#Separated part of MLP_p2 and MLP_p3
		self.norm_joint_scale = nn.LayerNorm(h_list_main[0])
		self.norm_joint_shift = nn.LayerNorm(h_list_main[0])

		self.linear_joint_scale = nn.Linear(h_list_joint2[-1], h_list_main[0])
		self.linear_joint_shift = nn.Linear(h_list_joint2[-1], h_list_main[0])

		self.dropout_joint_scale = nn.Dropout(p=dropout)
		self.dropout_joint_shift = nn.Dropout(p=dropout)



		self.dropout_att = nn.Dropout(p=dropout)

	def mask_softmax(self, A, mask):
		A_max = torch.max(A, dim=1, keepdim=True)[0]
		A_exp = torch.exp(A - A_max)
		A_exp = A_exp * mask # this step masks
		A_softmax = A_exp / torch.sum(A_exp ,dim=1, keepdim=True)
		return A_softmax


	def scale_and_shift(self, x, scale, shift):
		return (x * (scale + 1.0) + shift)


	#Feedforward preocess of MLP_p1
	def forward_joint1(self, x_joint):
		#feature processing before adding layer num info
		for i in range(self.num_layers_joint1):
			if i == 0:
				h_joint = self.linears_joint1[i](x_joint)
				h_joint = self.norm_joint1[i](h_joint)
				h_joint = self.dropout_joint1[i](h_joint)
				h_joint = F.gelu(h_joint)
			else:
				h_joint = self.linears_joint1[i](h_joint)
				h_joint = self.norm_joint1[i](h_joint)
				h_joint = self.dropout_joint1[i](h_joint)
				h_joint = F.gelu(h_joint)

		return h_joint

	#Feedforward preocess of MLP_p2 and MLP_p3
	def forward_joint2(self, h_joint, layer_num):

		time_emb = self.layer_num_emb(layer_num).view(1, -1)
		
		time_emb_scale = self.linear_layer_num_emb_scale(time_emb)
		time_emb_scale = self.dropout_layer_num_emb_scale(time_emb_scale)
		time_emb_scale = F.gelu(time_emb_scale)

		time_emb_shift = self.linear_layer_num_emb_shift(time_emb)
		time_emb_shift = self.dropout_layer_num_emb_shift(time_emb_shift)
		time_emb_shift = F.gelu(time_emb_shift)

		h_joint = self.scale_and_shift(h_joint, time_emb_scale, time_emb_shift)

		for i in range(self.num_layers_joint2):
			h_joint = self.linears_joint2[i](h_joint)
			h_joint = self.norm_joint2[i](h_joint)
			h_joint = self.dropout_joint2[i](h_joint)
			h_joint = F.gelu(h_joint)

		h_scale = self.linear_joint_scale(h_joint)
		h_scale = self.norm_joint_scale(h_scale)
		h_scale = self.dropout_joint_scale(h_scale)
		h_scale = F.gelu(h_scale)

		h_shift = self.linear_joint_shift(h_joint)
		h_shift = self.norm_joint_shift(h_shift)
		h_shift = self.dropout_joint_shift(h_shift)
		h_shift = F.gelu(h_shift)

		return h_scale, h_shift

	def att(self, key, query, value1, value2):

		att_score = torch.matmul(query, key.permute(1, 0))
		att_score = self.dropout_att(att_score)

		att_score = F.softmax(att_score, dim = -1)
		
		atted_value1 = torch.matmul(att_score, value1)
		atted_value2 = torch.matmul(att_score, value2)
		
		return atted_value1, atted_value2


	#Feedforward preocess of MLP_main
	def forward_main(self, x_main, PI):
		has_PI, x_main_mem, x_joint_mem = PI

		h_joint = self.forward_joint1(x_joint_mem)

		for i in range(self.num_layers_main):
			if i == 0:
				h_main = self.linears_main[i](x_main)
				h_main = self.norm_main[i](h_main)
				h_main = self.dropout_main[i](h_main)
				h_main = F.gelu(h_main)
			elif i == self.num_layers_main - 1:
				h_main = self.forest(h_main)
				return h_main[:, 1].view(-1, 1)
			else:
				h_main = self.linears_main[i](h_main)
				h_main = self.norm_main[i](h_main)
				h_main = self.dropout_main[i](h_main)
				h_main = F.gelu(h_main)

			if has_PI == False:
				if i == 0:
					h_tmp = self.linears_main[i](x_main_mem)
					h_tmp = self.norm_main[i](h_tmp)
					h_tmp = self.dropout_main[i](h_tmp)
					h_tmp = F.gelu(h_tmp)
				else:
					h_tmp = self.linears_main[i](h_tmp)
					h_tmp = self.norm_main[i](h_tmp)
					h_tmp = self.dropout_main[i](h_tmp)
					h_tmp = F.gelu(h_tmp)

			if has_PI:
				h_scale, h_shift = self.forward_joint2(h_joint, torch.FloatTensor(np.array([i])).cuda())
				h_main = self.scale_and_shift(h_main, h_scale, h_shift)
			else:
				h_scale, h_shift = self.forward_joint2(h_joint, torch.FloatTensor(np.array([i])).cuda())

				h_query = h_main
				h_key = h_tmp

				h_scale, h_shift = self.att(h_key, h_query, h_scale, h_shift)
				h_main = self.scale_and_shift(h_main, h_scale, h_shift)


class BCE_LOSS(nn.Module):
	def __init__(self):
		super().__init__()
		
	def forward(self, x, y):
		return F.binary_cross_entropy_with_logits(x, y)





