#ifndef GRNVIEWER_HPP
#define GRNVIEWER_HPP
#include <QQmlPropertyMap>
#include <QList>
#include <QObject>
#include <QAbstractListModel>
#include <QUrl>
#include <QTimer>
#include <QThread>
#include <QDebug>
#include <unordered_map>
#include <memory>
#include "../external/grgen/common.h"

class QProteinList : public QAbstractListModel {
	using QProtein = QPair<QString, qreal>;
	Q_OBJECT
 public:
	enum QProteinListRoles { NameRole = Qt::UserRole + 1, ValRole };
	QProteinList(QObject* = nullptr) {}
	void addProtein(QString n, qreal v) {
		beginInsertRows(QModelIndex(), rowCount(), rowCount());
		m_proteins << QProtein(n, v);
		endInsertRows();
	}
	int rowCount(const QModelIndex& parent = QModelIndex()) const {
		Q_UNUSED(parent);
		return m_proteins.count();
	};
	void clear() { removeRows(0, rowCount(), QModelIndex()); }
	bool removeRows(int position, int rows, const QModelIndex&) {
		beginRemoveRows(QModelIndex(), position, position + rows - 1);
		for (int row = 0; row < rows; ++row) {
			m_proteins.removeAt(position);
		}
		endRemoveRows();
		return true;
	}

	QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const {
		if (index.row() < 0 || index.row() >= m_proteins.count()) return QVariant();
		const QProtein& p = m_proteins[index.row()];
		if (role == NameRole)
			return p.first;
		else if (role == ValRole)
			return p.second;
		return QVariant();
	}

 private:
	QList<QProtein> m_proteins;
	QHash<int, QByteArray> roleNames() const {
		QHash<int, QByteArray> roles;
		roles[NameRole] = "name";
		roles[ValRole] = "val";
		return roles;
	}
};

class GrnUpdater : public QObject {
	Q_OBJECT
 public:
	GrnUpdater(QObject* = nullptr) {}

signals:
	void needUpdate();

 public slots:
	virtual void update(){};
};

template <typename GRN> struct GRNViewer : public GrnUpdater {
	QProteinList inputs, reguls, outputs;
	void* selected = nullptr;
	GRN grn;

	template <typename R> void preLoad(R* renderer) {
		auto engine = renderer->getEngine();
		engine->rootContext()->setContextProperty("grnInputs", &inputs);
		engine->rootContext()->setContextProperty("grnReguls", &reguls);
		engine->rootContext()->setContextProperty("grnOutputs", &outputs);
		engine->load(QUrl::fromLocalFile("../src/viewer/grnview.qml"));
		connect(this, SIGNAL(needUpdate()), this, SLOT(update()));
	}
	template <typename R> void postDraw(R* renderer) {
		if (renderer->getSelectedCell() != selected) {
			selected = renderer->getSelectedCell();
			if (selected)
				std::cerr << renderer->getSelectedCell()->ctrl.grn.toJSON() << std::endl;
		}
		if (selected) {
			grn = renderer->getSelectedCell()->ctrl.grn;
			emit needUpdate();
		}
	}
	virtual void update() {
		if (selected) {
			auto names = grn.getProteinNames(ProteinType::input);
			inputs.clear();
			for (auto& p : names) {
				QString pq = QString::fromStdString(p);
				inputs.addProtein(pq, grn.getProteinConcentration(p, ProteinType::input));
			}
			names = grn.getProteinNames(ProteinType::regul);
			reguls.clear();
			for (auto& p : names) {
				QString pq = QString::fromStdString(p);
				reguls.addProtein(pq, grn.getProteinConcentration(p, ProteinType::regul));
			}
			names = grn.getProteinNames(ProteinType::output);
			outputs.clear();
			for (auto& p : names) {
				QString pq = QString::fromStdString(p);
				outputs.addProtein(pq, grn.getProteinConcentration(p, ProteinType::output));
			}
		} else {
			inputs.clear();
			reguls.clear();
			outputs.clear();
		}
	}
};
#endif
