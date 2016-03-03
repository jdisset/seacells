#ifndef DRAWEXTENSIONS_HPP
#define DRAWEXTENSIONS_HPP
#include <mecacell/viewer/viewer.h>
#include <mecacell/viewer/primitives/deformablesphere.hpp>
#include <mecacell/viewer/primitives/deformabledisk.hpp>

using namespace MecacellViewer;
class ContactSurfaces {
	QOpenGLShaderProgram shader;
	DeformableDisk disk;

 public:
	ContactSurfaces() : disk(20) {}

	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/mvp.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/flat.frag"));
		shader.link();
		disk.load(shader);
	}

	template <typename R> void call(R *r) {
		using V = decltype(declval<typename R::Cell>().getPosition());
		const auto &view = r->getViewMatrix();
		const auto &projection = r->getProjectionMatrix();
		shader.bind();
		disk.vao.bind();
		QVector4D color(0.9f, 0.9f, 0.05f, 1.0f);
		shader.setUniformValue(shader.uniformLocation("color"), color);
		shader.setUniformValue(shader.uniformLocation("projection"), projection);
		shader.setUniformValue(shader.uniformLocation("view"), view);
		for (auto &con : r->getScenario().getWorld().cellCellConnections) {
			auto &cc = con.second;
			QMatrix4x4 model;
			model.translate(
			    toQV3D(cc.cells.first->getPosition() + cc.normal * cc.midpoint.first));
			auto rot = V::getRotation(V(0, 0, 1), cc.normal);
			model.rotate(rot.teta * 180.0 / M_PI, toQV3D(rot.n));
			float rad = static_cast<float>(sqrt(cc.sqradius));
			model.scale(rad, rad, rad);
			QMatrix4x4 nmatrix = (model).inverted().transposed();
			shader.setUniformValue(shader.uniformLocation("model"), model);
			shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
			GL->glDrawElements(GL_TRIANGLES, disk.indices.size(), GL_UNSIGNED_INT, 0);
		}
		color = QVector4D(0.1f, 0.7f, 0.f, 1.0f);
		shader.setUniformValue(shader.uniformLocation("color"), color);
		for (auto &con : r->getScenario().getWorld().cellCellConnections) {
			auto &cc = con.second;
			QMatrix4x4 model;
			model.translate(toQV3D(cc.cells.first->getPosition() +
			                       cc.icb.first.currentBasis.X * cc.midpoint.first));
			auto rot = V::getRotation(V(0, 0, 1), cc.icb.first.currentBasis.X);
			model.rotate(rot.teta * 180.0 / M_PI, toQV3D(rot.n));
			float rad = static_cast<float>(sqrt(cc.sqradius));
			model.scale(rad, rad, rad);
			QMatrix4x4 nmatrix = (model).inverted().transposed();
			shader.setUniformValue(shader.uniformLocation("model"), model);
			shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
			GL->glDrawElements(GL_TRIANGLES, disk.indices.size(), GL_UNSIGNED_INT, 0);
		}
		color = QVector4D(0.f, 0.1f, 0.7f, 1.0f);
		shader.setUniformValue(shader.uniformLocation("color"), color);
		for (auto &con : r->getScenario().getWorld().cellCellConnections) {
			auto &cc = con.second;
			QMatrix4x4 model;
			model.translate(toQV3D(cc.cells.second->getPosition() +
			                       cc.icb.second.currentBasis.X * cc.midpoint.second));
			auto rot = V::getRotation(V(0, 0, 1), cc.icb.second.currentBasis.X);
			model.rotate(rot.teta * 180.0 / M_PI, toQV3D(rot.n));
			float rad = static_cast<float>(sqrt(cc.sqradius));
			model.scale(rad, rad, rad);
			QMatrix4x4 nmatrix = (model).inverted().transposed();
			shader.setUniformValue(shader.uniformLocation("model"), model);
			shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
			GL->glDrawElements(GL_TRIANGLES, disk.indices.size(), GL_UNSIGNED_INT, 0);
		}
		disk.vao.release();
		shader.release();
	}
};
class ActiveConnections {
	QOpenGLShaderProgram shader;
	Cube cube;

 public:
	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/mvp.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/flat.frag"));
		shader.link();
		cube.load(shader);
	}

	template <typename R> void call(R *r) {
		const auto &view = r->getViewMatrix();
		const auto &projection = r->getProjectionMatrix();
		for (auto &con : r->getScenario().getWorld().cellCellConnections) {
			auto &cc = *(con.second.get());
			if (cc.adhCoef > 0) {
				shader.bind();
				cube.vao.bind();
				QColor color = QColor::fromHsvF(cc.adhCoef * 0.45, 0.7, 0.7);
				shader.setUniformValue(shader.uniformLocation("color"), color);
				shader.setUniformValue(shader.uniformLocation("projection"), projection);
				shader.setUniformValue(shader.uniformLocation("view"), view);
				QMatrix4x4 model;
				auto ab = toQV3D(cc.cells.second->getPosition() - cc.cells.first->getPosition());
				model.translate(toQV3D(cc.cells.second->getPosition()) - ab * 0.5);
				auto dp = ab.normalized().x();
				if (dp != 1 && dp != -1) {
					model.rotate(acos(dp) * 180.0 / M_PI,
					             QVector3D::crossProduct(QVector3D(1, 0, 0), ab));
					model.scale(ab.length() * 0.5, 1.0, 1.0);
					QMatrix4x4 nmatrix = (model).inverted().transposed();
					shader.setUniformValue(shader.uniformLocation("model"), model);
					shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
					GL->glDrawElements(GL_TRIANGLES, cube.indices.size(), GL_UNSIGNED_INT, 0);
				}
				cube.vao.release();
				shader.release();
			}
		}
	}
};
class DetailedConnections {
	QOpenGLShaderProgram shader;
	Cube cube;

 public:
	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/mvp.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/flat.frag"));
		shader.link();
		cube.load(shader);
	}

	template <typename R> void call(R *r) {
		const auto &view = r->getViewMatrix();
		const auto &projection = r->getProjectionMatrix();
		for (auto &con : r->getScenario().getWorld().cellCellConnections) {
			auto &cc = con.second;
			shader.bind();
			cube.vao.bind();
			QColor color = QColor::fromHsvF(0.7, 0.7, 0.7);
			shader.setUniformValue(shader.uniformLocation("color"), color);
			shader.setUniformValue(shader.uniformLocation("projection"), projection);
			shader.setUniformValue(shader.uniformLocation("view"), view);
			// first
			QMatrix4x4 model;
			auto ab =
			    toQV3D(cc.targets.first.b.X.rotated(cc.cells.first->getOrientationRotation()) *
			           cc.targets.first.d);
			model.translate(toQV3D(cc.cells.first->getPosition()) + ab * 0.5);
			auto dp = ab.normalized().x();
			if (dp != 1 && dp != -1) {
				model.rotate(acos(dp) * 180.0 / M_PI,
				             QVector3D::crossProduct(QVector3D(1, 0, 0), ab));
				model.scale(ab.length() * 0.5, 1.0, 1.0);
				QMatrix4x4 nmatrix = (model).inverted().transposed();
				shader.setUniformValue(shader.uniformLocation("model"), model);
				shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
				GL->glDrawElements(GL_TRIANGLES, cube.indices.size(), GL_UNSIGNED_INT, 0);
			}
			// second
			color = QColor::fromHsvF(0.3, 0.7, 0.7);
			shader.setUniformValue(shader.uniformLocation("color"), color);
			model = QMatrix4x4();
			ab = toQV3D(
			    cc.targets.second.b.X.rotated(cc.cells.second->getOrientationRotation()) *
			    cc.targets.second.d);
			model.translate(toQV3D(cc.cells.second->getPosition()) + ab * 0.5);
			dp = ab.normalized().x();
			if (dp != 1 && dp != -1) {
				model.rotate(acos(dp) * 180.0 / M_PI,
				             QVector3D::crossProduct(QVector3D(1, 0, 0), ab));
				model.scale(ab.length() * 0.5, 1.0, 1.0);
				QMatrix4x4 nmatrix = (model).inverted().transposed();
				shader.setUniformValue(shader.uniformLocation("model"), model);
				shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
				GL->glDrawElements(GL_TRIANGLES, cube.indices.size(), GL_UNSIGNED_INT, 0);
			}
			cube.vao.release();
			shader.release();
		}
	}
};

class Ground {
	QOpenGLShaderProgram shader;
	unique_ptr<QOpenGLTexture> texture = nullptr;
	Quad quad;

 public:
	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/texture.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/dumb.frag"));
		texture = unique_ptr<QOpenGLTexture>(
		    new QOpenGLTexture(QImage("../resources/sand.jpg").mirrored()));
		texture->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
		texture->setMagnificationFilter(QOpenGLTexture::Linear);
		shader.link();
		quad.load(shader);
	}

	template <typename R> void call(R *r) {
		const auto &view = r->getViewMatrix();
		const auto &projection = r->getProjectionMatrix();
		const float l = 100000.0f;
		shader.bind();
		quad.vao.bind();
		texture->bind(0);
		GL->glActiveTexture(GL_TEXTURE0);
		GL->glBindTexture(GL_TEXTURE_2D, texture->textureId());
		shader.setUniformValue(shader.uniformLocation("tex"), 0);
		shader.setUniformValue(shader.uniformLocation("projection"), projection);
		shader.setUniformValue(shader.uniformLocation("view"), view);
		shader.setUniformValue(shader.uniformLocation("texrepeat"), 300.0f);
		shader.setUniformValue(shader.uniformLocation("alpha"), 0.8f);
		QMatrix4x4 model;
		model.rotate(90.0f, QVector3D(1, 0, 0));
		model.translate(QVector3D(0, 0, 0));
		model.scale(l, l, l);
		shader.setUniformValue(shader.uniformLocation("model"), model);
		GL->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
		quad.vao.release();
		shader.release();
	}
};

class Nutrients {
	QOpenGLShaderProgram shader;
	unique_ptr<QOpenGLTexture> texture = nullptr;
	IcoSphere sphere;

 public:
	Nutrients() : sphere(1) {}
	void load() {
		shader.addShaderFromSourceCode(QOpenGLShader::Vertex,
		                               shaderWithHeader(":/shaders/mvp.vert"));
		shader.addShaderFromSourceCode(QOpenGLShader::Fragment,
		                               shaderWithHeader(":/shaders/flat.frag"));
		texture = unique_ptr<QOpenGLTexture>(
		    new QOpenGLTexture(QImage("../resources/sand.jpg").mirrored()));
		texture->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
		texture->setMagnificationFilter(QOpenGLTexture::Linear);
		shader.link();
		sphere.load(shader);
	}

	template <typename R> void call(R *r) {
		const auto &view = r->getViewMatrix();
		const auto &projection = r->getProjectionMatrix();
		shader.bind();
		sphere.vao.bind();
		texture->bind(0);
		shader.setUniformValue(shader.uniformLocation("projection"), projection);
		shader.setUniformValue(shader.uniformLocation("view"), view);
		for (auto &n : r->getScenario().nutrientSources) {
			QMatrix4x4 model;
			model.translate(n.pos.x(), n.pos.y(), n.pos.z());
			double c = n.content / n.initialcontent;
			double l = 15.0 + sqrt(n.sqradius * c) * 0.05;
			model.scale(l, l, l);
			QMatrix4x4 nmatrix = (model).inverted().transposed();
			shader.setUniformValue(shader.uniformLocation("model"), model);
			shader.setUniformValue(shader.uniformLocation("normalMatrix"), nmatrix);
			auto hsv = QColor::fromHsvF(c * 0.35, 0.9, 0.9);
			QVector4D col(hsv.redF(), hsv.greenF(), hsv.blueF(), 0.5);
			std::cerr << "c = " << c << ", r = " << col.x() << std::endl;
			shader.setUniformValue(shader.uniformLocation("color"), col);
			GL->glDrawElements(GL_TRIANGLES, sphere.indices.size(), GL_UNSIGNED_INT, 0);
		}
		sphere.vao.release();
		shader.release();
	}
};

#endif
